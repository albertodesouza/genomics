#!/usr/bin/env bash
#set -euo pipefail

#############################
# CONFIGURAÇÃO (AJUSTE AQUI)
#############################
# Genoma de referência (FASTA) e anotação (GTF)
REF=/dados/GENOMICS_DATA/top3/refs/reference.fa
GTF=/dados/GENOMICS_DATA/top3/refs/genes.gtf

# FASTQ entrada (R1). Aceita .fastq ou .fastq.gz
R1=/dados/GENOMICS_DATA/top3/fastq_ds/ERR3239334_1.ds.fastq.gz

# Threads e prefixo de saída
THREADS=128
OUT_PREFIX=ERR3239334_R1

# Strandedness para RNA-seq (bedtools). Opções: 0=desligado, 1=mesmo sentido, 2=oposto
# Para DNA-seq deixe STRANDED=0.
STRANDED=0

#############################
# CHECAGENS
#############################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERRO: '$1' não encontrado no PATH." >&2; exit 1; }; }
need bwa; need samtools; need bedtools; need awk; need sort; need grep; need cut; need comm

[[ -f "$REF" ]] || { echo "ERRO: REF não existe: $REF" >&2; exit 1; }
[[ -f "$GTF" ]] || { echo "ERRO: GTF não existe: $GTF" >&2; exit 1; }
[[ -f "$R1" ]] || { echo "ERRO: R1 não existe: $R1" >&2; exit 1; }

# Detecta R2 (paired) automaticamente
detect_r2() {
  local cand1 cand2
  cand1="${R1/_1/_2}"
  cand2="${R1/1.ds/2.ds}"
  if [[ -f "$cand1" ]]; then echo "$cand1"; return 0; fi
  if [[ -f "$cand2" ]]; then echo "$cand2"; return 0; fi
  echo ""  # não encontrou
}
R2="$(detect_r2)"
if [[ -n "$R2" ]]; then
  echo "[INFO] Detectado paired-end: R1=$R1  R2=$R2"
else
  echo "[INFO] Modo single-end: R1=$R1"
fi

# Lida com FASTQ.gz transparente
is_gz() { [[ "${1##*.}" == "gz" ]]; }

#############################
# 1) INDEXAR REFERÊNCIA (uma vez)
#############################
if [[ ! -f "${REF}.bwt" ]]; then
  echo "[INFO] Indexando referência para BWA..."
  bwa index "$REF"
fi
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

#############################
# 2) ALINHAMENTO + SORT + INDEX
#############################
BAM="${OUT_PREFIX}.sorted.bam"
if [[ -f "$BAM" ]]; then
  echo "[INFO] BAM já existe, pulando alinhamento: $BAM"
else
  echo "[INFO] Rodando alinhamento BWA-MEM..."
  if [[ -n "$R2" ]]; then
    # paired-end
    if is_gz "$R1"; then
      bwa mem -t "$THREADS" "$REF" <(gzip -cd "$R1") <(gzip -cd "$R2") \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$BAM"
    else
      bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$BAM"
    fi
  else
    # single-end
    if is_gz "$R1"; then
      bwa mem -t "$THREADS" "$REF" <(gzip -cd "$R1") \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$BAM"
    else
      bwa mem -t "$THREADS" "$REF" "$R1" \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$BAM"
    fi
  fi
fi
[[ -f "${BAM}.bai" ]] || samtools index "$BAM"

#############################
# 3) CONTAGEM POR CROMOSSOMO
#############################
CHROM_TSV="${OUT_PREFIX}.chrom_counts.tsv"
echo "[INFO] Gerando contagem por cromossomo -> $CHROM_TSV"
samtools idxstats "$BAM" | cut -f1,3 | grep -v '\*' > "$CHROM_TSV"
# colunas: chrom \t mapped_reads

#############################
# 4) GERAR BED DE GENES A PARTIR DO GTF
#############################
GENES_BED="genes.bed"
if [[ -f "$GENES_BED" ]]; then
  echo "[INFO] Usando genes.bed existente."
else
  echo "[INFO] Construindo genes.bed a partir de $GTF ..."
  awk 'BEGIN{FS=OFS="\t"}
  $3=="gene"{
    gid="."; gname=".";
    if (match($9,/gene_id "([^"]+)"/,a))   gid=a[1];
    if (match($9,/gene_name "([^"]+)"/,b)) gname=b[1];
    print $1,$4-1,$5,gid,gname,$7
  }' "$GTF" > "$GENES_BED"
fi
# Formato: chr  start  end  gene_id  gene_name  strand

#############################
# 5) READ -> (CHR, GENE) (MELHOR OVERLAP)
#############################
RAW_OV="${OUT_PREFIX}.read_gene_overlap.raw.bed"
BEST_TSV="${OUT_PREFIX}.read2gene.best.tsv"

echo "[INFO] BAM -> BED e interseção com genes..."
if [[ "$STRANDED" -eq 0 ]]; then
  bedtools bamtobed -i "$BAM" \
  | bedtools intersect -a - -b "$GENES_BED" -wo \
  > "$RAW_OV"
else
  # -s exige que o sentido da leitura corresponda ao gene (para RNA-seq com protocolo stranded)
  bedtools bamtobed -i "$BAM" \
  | bedtools intersect -a - -b "$GENES_BED" -wo -s \
  > "$RAW_OV"
fi

echo "[INFO] Selecionando gene de maior sobreposição por leitura..."
# Colunas RAW_OV:
# a.chr a.start a.end a.readName a.score a.strand  b.chr b.start b.end gene_id gene_name gene_strand overlap_bp
# Mantemos a primeira ocorrência após ordenar por readName e overlap desc.
sort -k4,4 -k13,13nr "$RAW_OV" \
| awk 'BEGIN{OFS="\t"} !seen[$4]++ { print $4,$1,$2,$3,$6,$10,$11,$13 }' \
> "$BEST_TSV"
# Saída: read_id chr start end read_strand gene_id gene_name overlap_bp

#############################
# 6) OPCIONAL: LISTAR INTERGÊNICAS E CONTAR POR GENE
#############################
ALL_IDS="${OUT_PREFIX}.all_mapped.ids"
GENE_IDS="${OUT_PREFIX}.in_genes.ids"
INTERGENIC_IDS="${OUT_PREFIX}.intergenic.ids"
GENE_COUNTS="${OUT_PREFIX}.gene_counts.tsv"

echo "[INFO] Listando leituras mapeadas e intergênicas..."
samtools view "$BAM" | cut -f1 | sort -u > "$ALL_IDS"
cut -f1 "$BEST_TSV" | sort -u > "$GENE_IDS"
comm -23 "$ALL_IDS" "$GENE_IDS" > "$INTERGENIC_IDS"

echo "[INFO] Contagem por gene -> $GENE_COUNTS"
cut -f6 "$BEST_TSV" | sort | uniq -c \
| awk 'BEGIN{OFS="\t"} {print $2,$1}' > "$GENE_COUNTS"
# colunas: gene_id \t n_reads

#############################
# RESUMO
#############################
echo "================= RESUMO ================="
echo "BAM ordenado:            $BAM"
echo "Index (BAI):             ${BAM}.bai"
echo "Contagem por cromossomo: $CHROM_TSV"
echo "Genes (BED):             $GENES_BED"
echo "Read→Gene (melhor):      $BEST_TSV"
echo "Intergênicas (ids):      $INTERGENIC_IDS"
echo "Contagem por gene:       $GENE_COUNTS"
echo "Stranded (0/1/2):        $STRANDED (aplicado no intersect)"
echo "Threads:                 $THREADS"
echo "=========================================="
