from __future__ import annotations

import math
from typing import Dict

import torch
import torch.nn as nn

from variant_transformer_predictor.config import ModelConfig
from variant_transformer_predictor.constants import BASE_PAD_ID
from variant_transformer_predictor.rope import apply_rope


class AlleleEncoder(nn.Module):
    def __init__(self, l_max: int, d_base: int, d_allele: int, dropout: float):
        super().__init__()
        self.l_max = l_max
        self.base_embedding = nn.Embedding(6, d_base, padding_idx=BASE_PAD_ID)
        self.network = nn.Sequential(
            nn.Linear(2 * l_max * d_base, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, d_allele),
        )

    def forward(self, ref_ids: torch.Tensor, alt_ids: torch.Tensor) -> torch.Tensor:
        ref = self.base_embedding(ref_ids).flatten(start_dim=-2)
        alt = self.base_embedding(alt_ids).flatten(start_dim=-2)
        return self.network(torch.cat([ref, alt], dim=-1))


class VariantTokenEncoder(nn.Module):
    def __init__(self, cfg: ModelConfig, num_genes: int, l_max: int):
        super().__init__()
        self.type_embedding = nn.Embedding(3, cfg.d_type)
        self.haplotype_embedding = nn.Embedding(2, cfg.d_hap)
        self.gene_embedding = nn.Embedding(num_genes, cfg.d_gene)
        self.length_projection = nn.Linear(1, cfg.d_len)
        self.allele_encoder = AlleleEncoder(l_max, cfg.d_base, cfg.d_allele, cfg.dropout)
        total = cfg.d_type + cfg.d_hap + cfg.d_gene + cfg.d_len + cfg.d_allele
        self.projection = nn.Linear(total, cfg.d_model)

    def forward(self, batch: Dict[str, torch.Tensor]) -> torch.Tensor:
        e_type = self.type_embedding(batch["variant_type"])
        e_hap = self.haplotype_embedding(batch["haplotype"])
        e_gene = self.gene_embedding(batch["gene"])
        e_len = self.length_projection(batch["length_norm"])
        e_allele = self.allele_encoder(batch["ref_allele"], batch["alt_allele"])
        return self.projection(torch.cat([e_type, e_hap, e_gene, e_len, e_allele], dim=-1))


class RoPEMultiheadSelfAttention(nn.Module):
    def __init__(self, d_model: int, heads: int, dropout: float, rope_base: float):
        super().__init__()
        self.d_model = d_model
        self.heads = heads
        self.head_dim = d_model // heads
        self.rope_base = rope_base
        self.qkv = nn.Linear(d_model, 3 * d_model)
        self.out = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, positions: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        bsz, seq_len, _ = x.shape
        qkv = self.qkv(x).view(bsz, seq_len, 3, self.heads, self.head_dim)
        q, k, v = qkv.unbind(dim=2)
        q = q.transpose(1, 2)
        k = k.transpose(1, 2)
        v = v.transpose(1, 2)
        q = apply_rope(q, positions, self.rope_base)
        k = apply_rope(k, positions, self.rope_base)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(self.head_dim)
        key_mask = attention_mask[:, None, None, :].to(dtype=torch.bool, device=x.device)
        scores = scores.masked_fill(~key_mask, torch.finfo(scores.dtype).min)
        attn = torch.softmax(scores, dim=-1)
        attn = self.dropout(attn)
        out = torch.matmul(attn, v).transpose(1, 2).contiguous().view(bsz, seq_len, self.d_model)
        return self.out(out)


class TransformerBlock(nn.Module):
    def __init__(self, cfg: ModelConfig):
        super().__init__()
        self.norm1 = nn.LayerNorm(cfg.d_model)
        self.attn = RoPEMultiheadSelfAttention(cfg.d_model, cfg.heads, cfg.dropout, cfg.rope_base)
        self.norm2 = nn.LayerNorm(cfg.d_model)
        hidden = cfg.d_model * cfg.mlp_ratio
        self.mlp = nn.Sequential(
            nn.Linear(cfg.d_model, hidden),
            nn.GELU(),
            nn.Dropout(cfg.dropout),
            nn.Linear(hidden, cfg.d_model),
            nn.Dropout(cfg.dropout),
        )

    def forward(self, x: torch.Tensor, positions: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        x = x + self.attn(self.norm1(x), positions, attention_mask)
        x = x + self.mlp(self.norm2(x))
        return x


class VariantTransformerClassifier(nn.Module):
    def __init__(self, cfg: ModelConfig, num_genes: int, num_classes: int, l_max: int):
        super().__init__()
        self.cfg = cfg
        self.token_encoder = VariantTokenEncoder(cfg, num_genes, l_max)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, cfg.d_model))
        self.dropout = nn.Dropout(cfg.dropout)
        self.blocks = nn.ModuleList([TransformerBlock(cfg) for _ in range(cfg.layers)])
        self.norm = nn.LayerNorm(cfg.d_model)
        self.classifier = nn.Sequential(
            nn.Linear(cfg.d_model, 256),
            nn.GELU(),
            nn.Dropout(cfg.dropout),
            nn.Linear(256, num_classes),
        )
        nn.init.normal_(self.cls_token, std=0.02)

    def forward(self, batch: Dict[str, torch.Tensor]) -> torch.Tensor:
        x = self.token_encoder(batch)
        bsz = x.shape[0]
        cls = self.cls_token.expand(bsz, -1, -1)
        x = torch.cat([cls, x], dim=1)
        x = self.dropout(x)
        zero = torch.zeros((bsz, 1), dtype=batch["position_relative"].dtype, device=x.device)
        positions = torch.cat([zero, batch["position_relative"]], dim=1)
        attention_mask = batch["attention_mask"].to(device=x.device, dtype=torch.bool)
        for block in self.blocks:
            x = block(x, positions, attention_mask)
        x = self.norm(x)
        return self.classifier(x[:, 0])
