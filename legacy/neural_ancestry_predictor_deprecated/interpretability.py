from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn
from rich.console import Console


console = Console()

# ==============================================================================
# INTERPRETABILITY: GRAD-CAM AND DEEPLIFT
# ==============================================================================

class GradCAM:
    """
    Grad-CAM (Gradient-weighted Class Activation Mapping) para CNNs.
    
    Calcula mapas de ativação que destacam regiões importantes para a predição
    de uma classe específica, usando gradientes da saída em relação às ativações
    da última camada convolucional.
    
    Referência: Selvaraju et al., 2017 - "Grad-CAM: Visual Explanations from 
    Deep Networks via Gradient-based Localization"
    
    Uso:
        gradcam = GradCAM(model, target_layer='auto')
        cam = gradcam.generate(input_tensor, target_class=None)  # None = classe predita
    """
    
    def __init__(self, model: nn.Module, target_layer: str = 'auto'):
        """
        Inicializa Grad-CAM.
        
        Args:
            model: Modelo CNN (CNNAncestryPredictor ou CNN2AncestryPredictor)
            target_layer: Camada alvo para extrair ativações
                         'auto' = detecta automaticamente a última conv
        """
        self.model = model
        self.target_layer = target_layer
        self.activations = None
        self.gradients = None
        self._hooks = []
        
        # Detectar tipo de modelo e camada alvo
        self._setup_target_layer()
    
    def _setup_target_layer(self):
        """Configura a camada alvo e registra hooks."""
        model_type = type(self.model).__name__
        
        if model_type == 'CNNAncestryPredictor':
            # CNN simples: usar conv1
            target = self.model.conv1
            self._layer_name = 'conv1'
        elif model_type == 'CNN2AncestryPredictor':
            # CNN2: usar última conv em features (índice 4 = Stage 3 Conv)
            # features = [Conv2d, ReLU, Conv2d, ReLU, Conv2d, ReLU]
            target = self.model.features[4]  # Stage 3 Conv2d
            self._layer_name = 'features[4] (Stage 3)'
        else:
            raise ValueError(f"Grad-CAM não suportado para modelo tipo: {model_type}")
        
        # Registrar hooks
        self._hooks.append(
            target.register_forward_hook(self._forward_hook)
        )
        self._hooks.append(
            target.register_full_backward_hook(self._backward_hook)
        )
        
        console.print(f"[cyan]Grad-CAM configurado para camada: {self._layer_name}[/cyan]")
    
    def _forward_hook(self, module, input, output):
        """Hook para capturar ativações no forward pass."""
        self.activations = output.detach()
    
    def _backward_hook(self, module, grad_input, grad_output):
        """Hook para capturar gradientes no backward pass."""
        self.gradients = grad_output[0].detach()
    
    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None) -> torch.Tensor:
        """
        Gera mapa Grad-CAM para uma entrada.
        
        Args:
            input_tensor: Tensor de entrada [batch, num_rows, effective_size]
            target_class: Classe alvo para gerar CAM. Se None, usa classe predita.
            
        Returns:
            Tensor com mapa de ativação [num_rows, effective_size] normalizado [0, 1]
        """
        # Garantir que modelo está em modo eval mas com gradientes
        was_training = self.model.training
        self.model.eval()
        
        # Forward pass
        input_tensor = input_tensor.requires_grad_(True)
        output = self.model(input_tensor)
        
        # Determinar classe alvo
        if target_class is None:
            target_class = output.argmax(dim=1).item()
        
        # Backward pass para a classe alvo
        # Estratégia inversa: propagar gradiente das outras classes para evitar
        # saturação quando a classe alvo tem probabilidade muito alta.
        # Gradiente = 0 para classe alvo, 1/(N-1) para as demais.
        # Isso captura "o que reduz as outras classes" ≈ "o que ativa a classe alvo"
        self.model.zero_grad()
        num_classes = output.shape[1]
        gradient = torch.ones_like(output) / (num_classes - 1)
        gradient[0, target_class] = 0
        output.backward(gradient=gradient, retain_graph=True)
        
        # Calcular pesos (média global dos gradientes por canal)
        # gradients shape: [batch, channels, height, width]
        weights = self.gradients.mean(dim=(2, 3), keepdim=True)  # [batch, channels, 1, 1]
        
        # Combinação ponderada das ativações
        # activations shape: [batch, channels, height, width]
        cam = (weights * self.activations).sum(dim=1, keepdim=True)  # [batch, 1, height, width]
        
        # ReLU para manter apenas contribuições positivas
        cam = torch.relu(cam)
        
        # Sem normalização - manter valores originais do CAM
        cam = cam.squeeze()  # [height, width]
        
        # Interpolar para tamanho da entrada original
        # Obter dimensão original da entrada (após seleção de genes)
        if hasattr(self.model, 'rows_to_use'):
            target_h = len(self.model.rows_to_use)
        else:
            target_h = input_tensor.shape[1]
        target_w = input_tensor.shape[2]
        
        # Redimensionar CAM para tamanho da entrada
        cam_resized = torch.nn.functional.interpolate(
            cam.unsqueeze(0).unsqueeze(0),
            size=(target_h, target_w),
            mode='bilinear',
            align_corners=False
        ).squeeze()
        
        # Restaurar estado do modelo
        if was_training:
            self.model.train()
        
        return cam_resized.cpu(), target_class
    
    def remove_hooks(self):
        """Remove todos os hooks registrados."""
        for hook in self._hooks:
            hook.remove()
        self._hooks = []


class DeepLIFT:
    """
    DeepLIFT (Deep Learning Important Features) para CNNs.
    
    Calcula atribuições de importância para cada feature de entrada,
    propagando diferenças de ativação (em relação a um baseline) através
    da rede usando regras de contribuição específicas.
    
    Implementação simplificada que usa gradientes × (input - baseline),
    uma aproximação válida para redes com ReLU (Rescale Rule).
    
    Referência: Shrikumar et al., 2017 - "Learning Important Features 
    Through Propagating Activation Differences"
    
    Uso:
        deeplift = DeepLIFT(model)
        attributions = deeplift.generate(input_tensor, baseline='zeros')
    """
    
    def __init__(self, model: nn.Module):
        """
        Inicializa DeepLIFT.
        
        Args:
            model: Modelo (NN, CNN ou CNN2)
        """
        self.model = model
        self._baseline_cache = None
        self._class_mean_cache: Dict[int, torch.Tensor] = {}  # Cache para média de atribuições por classe
        self._class_input_mean_cache: Dict[int, Tuple[torch.Tensor, int]] = {}  # Cache para (média das entradas, num_samples) por classe
    
    def _get_baseline(self, input_tensor: torch.Tensor, baseline_type: str,
                      dataset: Optional[Any] = None) -> torch.Tensor:
        """
        Obtém tensor baseline.
        
        Args:
            input_tensor: Tensor de entrada para obter shape e device
            baseline_type: 'zeros' ou 'mean'
            dataset: Dataset para calcular média (necessário se baseline_type='mean')
            
        Returns:
            Tensor baseline com mesmo shape que input_tensor
        """
        if baseline_type == 'zeros':
            return torch.zeros_like(input_tensor)
        
        elif baseline_type == 'mean':
            if self._baseline_cache is not None:
                return self._baseline_cache.to(input_tensor.device)
            
            if dataset is None:
                console.print("[yellow]⚠ Dataset não fornecido para baseline='mean', usando zeros[/yellow]")
                return torch.zeros_like(input_tensor)
            
            # Calcular média de todas as amostras
            console.print("\n\n[red]Calculando baseline (média do dataset)...[/red]\n\n")
            all_samples = []
            for i in range(min(len(dataset), 10000)):  # Limitar a 10000 amostras
                sample, _, _ = dataset[i]
                all_samples.append(sample)
            
            mean_baseline = torch.stack(all_samples).mean(dim=0)
            self._baseline_cache = mean_baseline
            return mean_baseline.unsqueeze(0).to(input_tensor.device)
        
        else:
            raise ValueError(f"Baseline type desconhecido: {baseline_type}")
    
    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None,
                 baseline_type: str = 'zeros', dataset: Optional[Any] = None) -> torch.Tensor:
        """
        Gera atribuições DeepLIFT para uma entrada.
        
        Usa a aproximação gradient × (input - baseline), que é equivalente
        ao DeepLIFT Rescale Rule para redes com ReLU.
        
        Args:
            input_tensor: Tensor de entrada [batch, num_rows, effective_size]
            target_class: Classe alvo. Se None, usa classe predita.
            baseline_type: 'zeros' ou 'mean'
            dataset: Dataset para calcular média (necessário se baseline_type='mean')
            
        Returns:
            Tensor de atribuições [num_rows, effective_size]
        """
        # Garantir modo eval
        was_training = self.model.training
        self.model.eval()
        
        # Obter baseline
        baseline = self._get_baseline(input_tensor, baseline_type, dataset)
        
        # Calcular diferença input - baseline
        delta = input_tensor - baseline
        
        # Habilitar gradientes
        input_tensor = input_tensor.clone().requires_grad_(True)
        
        # Forward pass
        output = self.model(input_tensor)
        
        # Determinar classe alvo
        if target_class is None:
            target_class = output.argmax(dim=1).item()
        
        # Backward pass
        self.model.zero_grad()
        one_hot = torch.zeros_like(output)
        one_hot[0, target_class] = 1
        output.backward(gradient=one_hot)
        
        # Atribuições = gradiente × delta (aproximação do DeepLIFT Rescale Rule)
        gradients = input_tensor.grad.detach()
        attributions = gradients * delta
        
        # Remover dimensão de batch
        attributions = attributions.squeeze(0)
        
        # Se modelo tem seleção de genes, aplicar máscara
        if hasattr(self.model, 'rows_to_use'):
            # Manter apenas as linhas usadas pelo modelo
            rows_to_use = self.model.rows_to_use
            # Criar tensor com zeros e preencher apenas as linhas usadas
            full_attr = torch.zeros_like(attributions)
            full_attr[rows_to_use, :] = attributions[rows_to_use, :]
            attributions = full_attr
        
        # Restaurar estado
        if was_training:
            self.model.train()
        
        return attributions.cpu(), target_class
    
    def generate_class_mean(self, target_class_idx: int, dataset: Any,
                            baseline_type: str = 'zeros') -> Tuple[torch.Tensor, torch.Tensor, int]:
        """
        Gera a média das atribuições DeepLIFT e das entradas para todas as amostras de uma classe.
        
        Args:
            target_class_idx: Índice da classe alvo
            dataset: Dataset para buscar amostras
            baseline_type: 'zeros' ou 'mean'
            
        Returns:
            Tupla (mean_attributions, mean_input, num_samples):
            - mean_attributions: Tensor com a média das atribuições [num_rows, effective_size]
            - mean_input: Tensor com a média das entradas [num_rows, effective_size]
            - num_samples: Número de amostras usadas
        """
        # Verificar cache
        if target_class_idx in self._class_mean_cache and target_class_idx in self._class_input_mean_cache:
            mean_input, num_samples = self._class_input_mean_cache[target_class_idx]
            return self._class_mean_cache[target_class_idx], mean_input, num_samples
        
        console.print(f"[cyan]Calculando média DeepLIFT para classe {target_class_idx}...[/cyan]")
        
        # Coletar todas as amostras da classe alvo
        class_samples = []
        for i in range(len(dataset)):
            sample, target, _ = dataset[i]
            if target == target_class_idx:
                class_samples.append(sample)
        
        if len(class_samples) == 0:
            console.print(f"[yellow]⚠ Nenhuma amostra encontrada para classe {target_class_idx}[/yellow]")
            return None, None, 0
        
        num_samples = len(class_samples)
        console.print(f"[cyan]  → {num_samples} amostras encontradas[/cyan]")
        
        # Calcular média das entradas
        mean_input = torch.stack(class_samples).mean(dim=0)
        
        # Calcular DeepLIFT para cada amostra
        all_attributions = []
        device = next(self.model.parameters()).device
        
        for i, sample in enumerate(class_samples):
            # Preparar entrada
            input_tensor = sample.unsqueeze(0).to(device)
            
            # Calcular atribuições para esta amostra
            attr, _ = self.generate(input_tensor, target_class=target_class_idx,
                                   baseline_type=baseline_type, dataset=dataset)
            all_attributions.append(attr)
            
            # Progress (a cada 10 amostras)
            if (i + 1) % 10 == 0:
                console.print(f"[dim]  → Processadas {i + 1}/{num_samples} amostras[/dim]")
        
        # Calcular média das atribuições
        mean_attributions = torch.stack(all_attributions).mean(dim=0)
        
        # Armazenar nos caches
        self._class_mean_cache[target_class_idx] = mean_attributions
        self._class_input_mean_cache[target_class_idx] = (mean_input, num_samples)
        
        console.print(f"[green]✓ Média DeepLIFT calculada para classe {target_class_idx}[/green]")
        
        return mean_attributions, mean_input, num_samples
    
    def find_max_individuals_for_regions(
        self, 
        top_regions: List[Tuple], 
        target_class_idx: int, 
        dataset: Any,
        baseline_type: str = 'zeros',
        tracks_per_gene: int = 6,
        gene_names: List[str] = None
    ) -> Dict[str, Dict]:
        """
        Encontra o indivíduo com maior valor DeepLIFT em cada região especificada.
        
        Args:
            top_regions: Lista de tuplas (gene_name, mean_val, chrom, genomic_pos, col_idx)
            target_class_idx: Índice da classe alvo
            dataset: Dataset para buscar amostras
            baseline_type: 'zeros' ou 'mean'
            tracks_per_gene: Número de tracks por gene (default 6)
            gene_names: Lista ordenada de nomes de genes
            
        Returns:
            Dict com estrutura:
            {
                gene_name: {
                    'sample_id': str,
                    'max_value': float,
                    'col_idx': int,
                    'all_region_values': Dict[str, float]  # valores em todas as 5 regiões
                }
            }
        """
        if not top_regions or gene_names is None:
            return {}
        
        console.print(f"\n[cyan]Buscando indivíduos com maior DeepLIFT em cada região...[/cyan]")
        
        # Coletar amostras da classe alvo
        class_samples = []
        class_sample_ids = []
        for i in range(len(dataset)):
            sample, target, sample_idx = dataset[i]
            if target == target_class_idx:
                # Obter sample_id e metadata completo
                if hasattr(dataset, 'get_sample_metadata'):
                    sample_meta = dataset.get_sample_metadata(i)
                    sample_id = sample_meta.get('sample_id', f'sample_{i}')
                    superpopulation = sample_meta.get('superpopulation', 'UNK')
                    population = sample_meta.get('population', 'UNK')
                elif hasattr(dataset, '_sample_metadata') and dataset._sample_metadata:
                    sample_id = dataset._sample_metadata.get(i, {}).get('sample_id', f'sample_{i}')
                    superpopulation = 'UNK'
                    population = 'UNK'
                else:
                    sample_id = f'sample_{i}'
                    superpopulation = 'UNK'
                    population = 'UNK'
                class_samples.append((sample, i, sample_id, superpopulation, population))
                class_sample_ids.append(sample_id)
        
        if not class_samples:
            console.print(f"[yellow]⚠ Nenhuma amostra encontrada para classe {target_class_idx}[/yellow]")
            return {}
        
        console.print(f"[cyan]  → {len(class_samples)} amostras da classe alvo[/cyan]")
        
        # Para cada região, encontrar o indivíduo com maior valor
        device = next(self.model.parameters()).device
        results = {}
        
        # Cache de atribuições individuais (sample_idx -> attributions)
        attribution_cache = {}
        
        num_regions = len(top_regions)
        for region_idx, region_tuple in enumerate(top_regions):
            # Suportar tanto formato per_gene (5 campos) quanto global (7 campos)
            gene_name = region_tuple[0]
            mean_val = region_tuple[1]
            chrom = region_tuple[2]
            genomic_pos = region_tuple[3]
            col_idx = region_tuple[4]
            track_idx = region_tuple[6] if len(region_tuple) > 6 else None
            
            # Chave única para modo global (diferencia múltiplas entradas do mesmo gene)
            region_key = f"{gene_name}_{col_idx}_{track_idx}" if track_idx is not None else gene_name
            
            track_info = f" (track {track_idx})" if track_idx is not None else ""
            console.print(f"[dim]  Processando região {region_idx + 1}/{num_regions}: {gene_name}{track_info}...[/dim]")
            
            # Encontrar índice do gene
            if gene_name not in gene_names:
                continue
            gene_idx = gene_names.index(gene_name)
            start_row = gene_idx * tracks_per_gene
            end_row = (gene_idx + 1) * tracks_per_gene
            
            max_value = -float('inf')
            max_sample_id = None
            max_sample_cache_idx = None
            max_superpopulation = 'UNK'
            max_population = 'UNK'
            
            for sample, sample_cache_idx, sample_id, superpopulation, population in class_samples:
                # Verificar cache
                if sample_cache_idx not in attribution_cache:
                    input_tensor = sample.unsqueeze(0).to(device)
                    attr, _ = self.generate(input_tensor, target_class=target_class_idx,
                                           baseline_type=baseline_type, dataset=dataset)
                    attribution_cache[sample_cache_idx] = attr.numpy()
                
                attr_np = attribution_cache[sample_cache_idx]
                
                # Extrair valor na região específica
                if end_row <= attr_np.shape[0] and col_idx < attr_np.shape[1]:
                    # Média dos valores nas tracks do gene na coluna específica
                    region_value = attr_np[start_row:end_row, col_idx].mean()
                    
                    if region_value > max_value:
                        max_value = region_value
                        max_sample_id = sample_id
                        max_sample_cache_idx = sample_cache_idx
                        max_superpopulation = superpopulation
                        max_population = population
            
            if max_sample_id is not None:
                # Calcular valores em todas as N regiões para este indivíduo
                all_region_values = {}
                attr_np = attribution_cache[max_sample_cache_idx]
                
                for other_region in top_regions:
                    other_gene = other_region[0]
                    other_col_idx = other_region[4]
                    other_track_idx = other_region[6] if len(other_region) > 6 else None
                    other_key = f"{other_gene}_{other_col_idx}_{other_track_idx}" if other_track_idx is not None else other_gene
                    
                    if other_gene in gene_names:
                        other_gene_idx = gene_names.index(other_gene)
                        other_start = other_gene_idx * tracks_per_gene
                        other_end = (other_gene_idx + 1) * tracks_per_gene
                        if other_end <= attr_np.shape[0] and other_col_idx < attr_np.shape[1]:
                            all_region_values[other_key] = float(attr_np[other_start:other_end, other_col_idx].mean())
                
                results[region_key] = {
                    'sample_id': max_sample_id,
                    'superpopulation': max_superpopulation,
                    'population': max_population,
                    'max_value': float(max_value),
                    'col_idx': col_idx,
                    'genomic_pos': genomic_pos,
                    'chrom': chrom,
                    'all_region_values': all_region_values,
                    'gene_name': gene_name,
                    'track_idx': track_idx
                }
                
                console.print(f"[green]    ✓ {gene_name}{track_info}: {max_sample_id} ({max_superpopulation}/{max_population}, valor = {max_value:.6f})[/green]")
        
        return results


