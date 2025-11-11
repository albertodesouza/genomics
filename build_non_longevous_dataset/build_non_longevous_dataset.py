#!/usr/bin/env python3
"""
build_non_longevous_dataset.py

Pipeline to build datasets from non-longevous individuals from the 1000 Genomes Project.

This program:
1. Analyzes a CSV file with metadata about non-longevous individuals
2. Prints statistics about superpopulations, populations, and sex distribution
3. Selects samples according to user-specified criteria
4. Runs build_window_and_predict.py for each selected sample

Requirements:
  - Python 3.8+
  - pandas
  - pyyaml
  - build_window_and_predict.py (in the parent directory)

CSV Format:
  Header: FamilyID,SampleID,FatherID,MotherID,Sex,Population,Superpopulation
  Sex: 1=Male, 2=Female
  Example:
    BB01,HG01879,0,0,1,ACB,AFR
    BB01,HG01880,0,0,2,ACB,AFR
    BB01,HG01881,HG01879,HG01880,2,ACB,AFR

Usage:
  # Step 1: Analyze metadata and print statistics
  python3 build_non_longevous_dataset.py --config configs/default.yaml
  
  # Step 2: Enable sample selection and run predictions
  # (Edit configs/default.yaml to enable additional steps)
  python3 build_non_longevous_dataset.py --config configs/default.yaml

Author: ChatGPT (for Alberto)
Last updated: 2025-11-04
"""

import argparse
import os
import sys
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict, Counter
import yaml
import pandas as pd
import numpy as np

# Import dos novos módulos para Dataset PyTorch
try:
    from .frog_ancestry_parser import FROGAncestryData
    from .dataset_builder import IndividualDatasetBuilder, DatasetMetadataBuilder
except ImportError:
    # Fallback para execução como script standalone
    from frog_ancestry_parser import FROGAncestryData
    from dataset_builder import IndividualDatasetBuilder, DatasetMetadataBuilder


def load_config(config_path: Path) -> dict:
    """
    Load YAML configuration file.
    
    Relative paths in the config file are resolved relative to the config file's directory.
    """
    config_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Expand environment variables and resolve relative paths
    def expand_env_vars(obj, is_path_field=False):
        if isinstance(obj, dict):
            return {k: expand_env_vars(v, k in ['metadata_csv', 'fasta', 'vcf_pattern']) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [expand_env_vars(item) for item in obj]
        elif isinstance(obj, str):
            # Expand environment variables
            if obj.startswith("${") and obj.endswith("}"):
                env_var = obj[2:-1]
                obj = os.environ.get(env_var, obj)
            
            # Resolve relative paths (only for known path fields)
            if is_path_field and not obj.startswith("/"):
                obj = str((config_dir / obj).resolve())
            
            return obj
        return obj
    
    return expand_env_vars(config)


def load_metadata_csv(csv_path: Path) -> pd.DataFrame:
    """Load 1000 Genomes metadata CSV."""
    print(f"[INFO] Carregando arquivo CSV: {csv_path}")
    
    if not csv_path.exists():
        raise FileNotFoundError(f"Arquivo CSV não encontrado: {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Validate required columns
    required_cols = ['FamilyID', 'SampleID', 'FatherID', 'MotherID', 
                     'Sex', 'Population', 'Superpopulation']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        raise ValueError(f"Colunas ausentes no CSV: {missing_cols}")
    
    print(f"[INFO] CSV carregado: {len(df)} indivíduos")
    return df


def analyze_metadata(df: pd.DataFrame) -> Dict:
    """
    Analyze metadata and return statistics.
    
    Returns:
        Dict with statistics about superpopulations, populations, and sex distribution
    """
    stats = {
        'total_samples': len(df),
        'superpopulations': {},
        'populations': {}
    }
    
    # Group by superpopulation
    for superpop in sorted(df['Superpopulation'].unique()):
        superpop_df = df[df['Superpopulation'] == superpop]
        
        # Count populations in this superpopulation
        populations_in_superpop = superpop_df['Population'].unique()
        
        # Count sex distribution in superpopulation
        sex_counts = superpop_df['Sex'].value_counts().to_dict()
        
        stats['superpopulations'][superpop] = {
            'total_samples': len(superpop_df),
            'n_populations': len(populations_in_superpop),
            'populations': list(populations_in_superpop),
            'sex_distribution': {
                'male': sex_counts.get(1, 0),
                'female': sex_counts.get(2, 0)
            }
        }
    
    # Group by population
    for pop in sorted(df['Population'].unique()):
        pop_df = df[df['Population'] == pop]
        superpop = pop_df['Superpopulation'].iloc[0]
        
        # Count sex distribution in population
        sex_counts = pop_df['Sex'].value_counts().to_dict()
        
        stats['populations'][pop] = {
            'superpopulation': superpop,
            'total_samples': len(pop_df),
            'sex_distribution': {
                'male': sex_counts.get(1, 0),
                'female': sex_counts.get(2, 0)
            }
        }
    
    return stats


def print_statistics(stats: Dict):
    """Print formatted statistics."""
    print("\n" + "="*80)
    print("ESTATÍSTICAS DO DATASET - 1000 GENOMES PROJECT")
    print("="*80)
    
    print(f"\n📊 TOTAL DE AMOSTRAS: {stats['total_samples']}")
    
    # Superpopulations
    print(f"\n🌍 SUPERPOPULAÇÕES: {len(stats['superpopulations'])}")
    print("-" * 80)
    
    for superpop, data in sorted(stats['superpopulations'].items()):
        print(f"\n  {superpop}:")
        print(f"    • Total de indivíduos: {data['total_samples']}")
        print(f"    • Masculino: {data['sex_distribution']['male']}")
        print(f"    • Feminino: {data['sex_distribution']['female']}")
        print(f"    • Número de populações: {data['n_populations']}")
        print(f"    • Populações: {', '.join(data['populations'])}")
    
    # Populations
    print(f"\n🏘️  POPULAÇÕES: {len(stats['populations'])}")
    print("-" * 80)
    
    # Group populations by superpopulation for better display
    pops_by_superpop = defaultdict(list)
    for pop, data in stats['populations'].items():
        pops_by_superpop[data['superpopulation']].append((pop, data))
    
    for superpop in sorted(pops_by_superpop.keys()):
        print(f"\n  {superpop}:")
        for pop, data in sorted(pops_by_superpop[superpop]):
            print(f"    {pop}: {data['total_samples']} indivíduos "
                  f"(♂ {data['sex_distribution']['male']}, "
                  f"♀ {data['sex_distribution']['female']})")
    
    print("\n" + "="*80)
    print("LEGENDA:")
    print("  ♂ = Masculino (Sex=1)")
    print("  ♀ = Feminino (Sex=2)")
    print("="*80 + "\n")


def select_samples(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """
    Select samples according to configuration criteria.
    
    Returns:
        DataFrame with selected samples
    """
    sel_config = config['sample_selection']
    
    level = sel_config['level']  # "superpopulation" or "population"
    samples_per_group = sel_config['samples_per_group']
    include_groups = sel_config.get('include_groups', [])
    exclude_groups = sel_config.get('exclude_groups', [])
    sex_filter = sel_config.get('sex_filter', 'all')
    random_seed = sel_config.get('random_seed', 42)
    
    print(f"\n[INFO] Selecionando amostras...")
    print(f"  • Nível: {level}")
    print(f"  • Amostras por grupo: {samples_per_group}")
    
    # Apply sex filter
    if sex_filter == 'male':
        df = df[df['Sex'] == 1]
        print(f"  • Filtro de sexo: apenas masculino")
    elif sex_filter == 'female':
        df = df[df['Sex'] == 2]
        print(f"  • Filtro de sexo: apenas feminino")
    else:
        print(f"  • Filtro de sexo: todos")
    
    # Determine grouping column
    group_col = 'Superpopulation' if level == 'superpopulation' else 'Population'
    
    # Filter groups
    if include_groups:
        df = df[df[group_col].isin(include_groups)]
        print(f"  • Incluindo apenas: {', '.join(include_groups)}")
    
    if exclude_groups:
        df = df[~df[group_col].isin(exclude_groups)]
        print(f"  • Excluindo: {', '.join(exclude_groups)}")
    
    # Sample from each group
    selected_samples = []
    
    for group in sorted(df[group_col].unique()):
        group_df = df[df[group_col] == group]
        
        # Sample up to samples_per_group
        n_to_sample = min(samples_per_group, len(group_df))
        sampled = group_df.sample(n=n_to_sample, random_state=random_seed)
        
        selected_samples.append(sampled)
        print(f"  • {group}: selecionados {len(sampled)} de {len(group_df)} disponíveis")
    
    result_df = pd.concat(selected_samples, ignore_index=True)
    
    print(f"\n[INFO] Total de amostras selecionadas: {len(result_df)}")
    
    return result_df


def load_checkpoint(checkpoint_file: Path) -> dict:
    """Load checkpoint file if it exists."""
    if checkpoint_file.exists():
        with open(checkpoint_file, 'r') as f:
            return json.load(f)
    return {'completed_samples': [], 'failed_samples': []}


def save_checkpoint(checkpoint_file: Path, checkpoint_data: dict):
    """Save checkpoint file."""
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoint_data, f, indent=2)


def get_chromosome_for_gene(config: dict) -> str:
    """
    Determine which chromosome contains the target gene.
    This is a simplified version - in production you'd query the GTF.
    For now, we'll assume it's specified or we detect it.
    """
    # This could be enhanced to actually parse GTF and find the chromosome
    # For now, return a placeholder that will be determined at runtime
    return None


def validate_vcf_exists(sample_id: str, vcf_pattern: str, chromosome: str) -> Optional[Path]:
    """
    Check if VCF file exists for the given sample and chromosome.
    
    Returns:
        Path to VCF if it exists and contains the sample, None otherwise
    """
    vcf_path = Path(vcf_pattern.format(chrom=chromosome))
    
    if not vcf_path.exists():
        print(f"[WARN] VCF não encontrado: {vcf_path}")
        return None
    
    # Check if VCF index exists
    vcf_idx = Path(str(vcf_path) + ".tbi")
    if not vcf_idx.exists():
        print(f"[WARN] Índice VCF não encontrado: {vcf_idx}")
        return None
    
    # Optionally, we could check if sample is in VCF using bcftools query -l
    # For now, assume it exists if file exists
    
    return vcf_path


def run_build_window_predict(
    sample_id: str, 
    config: dict, 
    output_dir: Path,
    target_name: Optional[str] = None
) -> Tuple[bool, Optional[str]]:
    """
    Run build_window_and_predict.py for a single sample.
    
    Args:
        sample_id: ID do sample
        config: Configuração do pipeline
        output_dir: Diretório de saída (será usado como base para individuals/)
        target_name: Nome do gene/SNP sendo processado (para logging)
    
    Returns:
        Tupla (success: bool, target_name: Optional[str])
        target_name é o nome do gene/SNP processado
    """
    params = config['build_window_params']
    data_sources = config['data_sources']
    
    # Build command
    # build_window_and_predict.py está no diretório pai
    script_dir = Path(__file__).parent
    bwp_script = script_dir.parent / "build_window_and_predict.py"
    
    # Determinar output directory: individuals/SAMPLE_ID/windows/
    # O build_window_and_predict.py criará SAMPLE__TARGET dentro deste dir
    # Então precisamos apontar para individuals/SAMPLE_ID/ e ele criará windows/
    individual_base = output_dir / "individuals" / sample_id
    individual_base.mkdir(parents=True, exist_ok=True)
    
    # O output será: individuals/SAMPLE_ID/SAMPLE__TARGET/
    # Precisamos mover depois para: individuals/SAMPLE_ID/windows/TARGET/
    
    cmd = [
        sys.executable,  # Use same Python interpreter
        str(bwp_script),
        "--sample", sample_id,
        "--ref-fasta", str(Path(data_sources['reference']['fasta']).expanduser()),
        "--outdir", str(individual_base),
    ]
    
    # Add mode
    mode = params.get('mode', 'gene')
    cmd.extend(["--mode", mode])
    
    # Mode-specific arguments
    if mode == 'gene':
        gene_config = params.get('gene', {})
        
        # Check for gene list file first
        if gene_config.get('gene_list_file'):
            cmd.extend(["--gene-list-file", gene_config['gene_list_file']])
        # Otherwise check for single gene
        elif gene_config.get('symbol'):
            cmd.extend(["--gene", gene_config['symbol']])
        elif gene_config.get('id'):
            cmd.extend(["--gene-id", gene_config['id']])
        else:
            raise ValueError("Gene mode requires either gene.symbol, gene.id, or gene.gene_list_file in config")
    
    elif mode == 'snp':
        snp_config = params.get('snp', {})
        
        if not snp_config.get('snp_list_file'):
            raise ValueError("SNP mode requires snp.snp_list_file in config")
        
        cmd.extend(["--snp-list-file", snp_config['snp_list_file']])
    
    # Add GTF if specified
    if params.get('gtf_feather'):
        cmd.extend(["--gtf-feather", params['gtf_feather']])
    
    # Add window size
    if params.get('window_size'):
        cmd.extend(["--window-size", str(params['window_size'])])
    
    # Add haplotype options
    if params.get('skip_h2'):
        cmd.append("--skip-h2")
    
    if params.get('also_iupac'):
        cmd.append("--also-iupac")
    
    # Add prediction options
    if params.get('predict'):
        cmd.append("--predict")
        
        # Add API key
        if params.get('api_key'):
            cmd.extend(["--api-key", params['api_key']])
        
        # Add outputs
        if params.get('outputs'):
            cmd.extend(["--outputs", params['outputs']])
        
        # Add ontology
        if params.get('ontology'):
            cmd.extend(["--ontology", params['ontology']])
        
        # Add all-tissues flag
        if params.get('all_tissues'):
            cmd.append("--all-tissues")
    
    # VCF path - we need to determine which chromosome
    # For simplicity, we'll add a placeholder and let build_window_and_predict handle it
    # In a real implementation, you'd need to know which chromosome contains your gene
    
    # IMPORTANT: We need to determine the chromosome for the VCF
    # This requires either:
    # 1. Loading GTF and finding the gene
    # 2. Having the user specify it in config
    # For now, we'll assume a single VCF or require chromosome in config
    
    # Check if there's a chromosome specified in config
    vcf_pattern = data_sources.get('vcf_pattern', '')
    
    # If pattern contains {chrom}, we need to determine chromosome
    if '{chrom}' in vcf_pattern:
        print(f"[WARN] VCF pattern contém {{chrom}}, mas cromossomo não foi determinado.")
        print(f"[WARN] Por favor, especifique o cromossomo no config ou forneça VCF completo.")
        # For now, try common chromosomes or skip
        # This is a limitation that should be addressed in production
        return False
    else:
        cmd.extend(["--vcf", vcf_pattern])
    
    # Determinar target_name para retorno
    extracted_target_name = None
    mode = params.get('mode', 'gene')
    
    if mode == 'gene':
        gene_config = params.get('gene', {})
        if gene_config.get('symbol'):
            extracted_target_name = gene_config['symbol']
        elif gene_config.get('id'):
            extracted_target_name = gene_config['id']
    elif mode == 'snp':
        # Para SNP mode, será criado múltiplos diretórios (um por SNP)
        # Retornaremos None e o caller precisa descobrir quais foram criados
        extracted_target_name = None
    
    # Run command
    print(f"\n[INFO] Executando: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print(result.stdout)
        
        # Reorganizar estrutura de saída
        # build_window_and_predict.py cria: outdir/SAMPLE__TARGET/
        # Queremos mover para: individuals/SAMPLE/windows/TARGET/
        _reorganize_output_structure(individual_base, sample_id)
        
        return True, extracted_target_name
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Falha ao processar amostra {sample_id}")
        print(f"[ERROR] Código de saída: {e.returncode}")
        print(f"[ERROR] stderr: {e.stderr}")
        return False, None


def _reorganize_output_structure(individual_base: Path, sample_id: str) -> List[str]:
    """
    Reorganiza estrutura de saída do build_window_and_predict.py.
    
    Move de: individual_base/SAMPLE__TARGET/
    Para: individual_base/windows/TARGET/
    
    Args:
        individual_base: Diretório base do indivíduo (individuals/SAMPLE_ID/)
        sample_id: ID do sample
        
    Returns:
        Lista de nomes de targets processados
    """
    windows_dir = individual_base / "windows"
    windows_dir.mkdir(exist_ok=True)
    
    # Encontrar diretórios SAMPLE__TARGET
    pattern = f"{sample_id}__*"
    target_dirs = list(individual_base.glob(pattern))
    
    processed_targets = []
    
    for target_dir in target_dirs:
        if not target_dir.is_dir():
            continue
        
        # Extrair nome do target
        target_name = target_dir.name.replace(f"{sample_id}__", "")
        
        # Destino final
        final_dir = windows_dir / target_name
        
        # Mover (ou sobrescrever se já existir)
        if final_dir.exists():
            import shutil
            shutil.rmtree(final_dir)
        
        target_dir.rename(final_dir)
        processed_targets.append(target_name)
        
        print(f"[INFO] Janela organizada: {target_name}")
    
    return processed_targets


def generate_summary_report(config: dict, checkpoint_data: dict, output_dir: Path):
    """Generate a summary report of the processing."""
    report_file = output_dir / "processing_summary.txt"
    
    completed = checkpoint_data.get('completed_samples', [])
    failed = checkpoint_data.get('failed_samples', [])
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("NON-LONGEVOUS DATASET - RELATÓRIO DE PROCESSAMENTO\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Total de amostras completadas: {len(completed)}\n")
        f.write(f"Total de amostras com falha: {len(failed)}\n\n")
        
        if completed:
            f.write("Amostras completadas:\n")
            for sample in completed:
                f.write(f"  • {sample}\n")
            f.write("\n")
        
        if failed:
            f.write("Amostras com falha:\n")
            for sample in failed:
                f.write(f"  • {sample}\n")
            f.write("\n")
        
        f.write("="*80 + "\n")
    
    print(f"\n[INFO] Relatório salvo em: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Build non-longevous dataset from 1000 Genomes Project"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="configs/default.yaml",
        help="Path to YAML configuration file (default: configs/default.yaml)"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"[ERROR] Arquivo de configuração não encontrado: {config_path}")
        sys.exit(1)
    
    config = load_config(config_path)
    print(f"[INFO] Configuração carregada: {config_path}")
    
    # Create output directory
    output_dir = Path(config['project']['output_dir']).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Diretório de saída: {output_dir}")
    
    # Get pipeline steps
    steps = config['pipeline']['steps']
    
    # Step 1: Analyze metadata
    if steps.get('analyze_metadata', False):
        print("\n" + "="*80)
        print("PASSO 1: ANÁLISE DE METADADOS")
        print("="*80)
        
        csv_path = Path(config['data_sources']['metadata_csv']).resolve()
        df = load_metadata_csv(csv_path)
        
        stats = analyze_metadata(df)
        print_statistics(stats)
        
        # Save statistics to JSON
        stats_file = output_dir / "metadata_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"[INFO] Estatísticas salvas em: {stats_file}")
    
    # Step 2: Select samples
    selected_df = None
    if steps.get('select_samples', False):
        print("\n" + "="*80)
        print("PASSO 2: SELEÇÃO DE AMOSTRAS")
        print("="*80)
        
        if df is None:
            csv_path = Path(config['data_sources']['metadata_csv']).resolve()
            df = load_metadata_csv(csv_path)
        
        selected_df = select_samples(df, config)
        
        # Save selected samples
        selected_file = output_dir / "selected_samples.csv"
        selected_df.to_csv(selected_file, index=False)
        print(f"[INFO] Amostras selecionadas salvas em: {selected_file}")
    
    # Step 3: Validate VCFs (optional)
    if steps.get('validate_vcfs', False):
        print("\n" + "="*80)
        print("PASSO 3: VALIDAÇÃO DE ARQUIVOS VCF")
        print("="*80)
        
        if selected_df is None:
            selected_file = output_dir / "selected_samples.csv"
            if selected_file.exists():
                selected_df = pd.read_csv(selected_file)
            else:
                print("[ERROR] Nenhuma amostra selecionada. Execute o passo 'select_samples' primeiro.")
                sys.exit(1)
        
        # This step would validate VCF existence
        # Skipped for now due to chromosome determination complexity
        print("[INFO] Validação de VCF não implementada nesta versão.")
        print("[INFO] VCFs serão validados durante a execução de build_window_and_predict.py")
    
    # Step 4: Run predictions
    if steps.get('run_predictions', False):
        print("\n" + "="*80)
        print("PASSO 4: EXECUÇÃO DE PREDIÇÕES")
        print("="*80)
        
        if selected_df is None:
            selected_file = output_dir / "selected_samples.csv"
            if selected_file.exists():
                selected_df = pd.read_csv(selected_file)
            else:
                print("[ERROR] Nenhuma amostra selecionada. Execute o passo 'select_samples' primeiro.")
                sys.exit(1)
        
        # Carregar FROGAncestryCalc likelihoods
        print("\n[INFO] Carregando dados de ancestralidade do FROGAncestryCalc...")
        frog_data = None
        try:
            script_dir = Path(__file__).parent
            project_root = script_dir.parent
            likelihood_file = project_root / "FROGAncestryCalc/output/whole_1000genomes_55aisnps_likelihood.txt"
            mapping_file = project_root / "FROGAncestryCalc/population_mapping_1000genomes.csv"
            
            if likelihood_file.exists() and mapping_file.exists():
                frog_data = FROGAncestryData(likelihood_file, mapping_file)
            else:
                print("[WARN] Arquivos do FROGAncestryCalc não encontrados. Continuando sem dados de ancestralidade.")
        except Exception as e:
            print(f"[WARN] Erro ao carregar FROGAncestryCalc: {e}")
            print("[WARN] Continuando sem dados de ancestralidade.")
        
        # Load checkpoint
        checkpoint_file = output_dir / config['pipeline']['checkpoint_file']
        checkpoint_data = load_checkpoint(checkpoint_file)
        
        # Obter parâmetros de build_window
        params = config['build_window_params']
        
        # Process each sample
        total_samples = len(selected_df)
        for idx, row in selected_df.iterrows():
            sample_id = row['SampleID']
            
            # Skip if already completed
            if sample_id in checkpoint_data['completed_samples']:
                print(f"\n[INFO] Amostra {sample_id} já processada (passo {idx+1}/{total_samples})")
                continue
            
            print(f"\n[INFO] Processando amostra {sample_id} ({idx+1}/{total_samples})")
            print(f"  • População: {row['Population']}")
            print(f"  • Superpopulação: {row['Superpopulation']}")
            print(f"  • Sexo: {'Masculino' if row['Sex'] == 1 else 'Feminino'}")
            
            # Preparar informações do sample
            sample_info = {
                'FamilyID': row.get('FamilyID', '0'),
                'SampleID': sample_id,
                'Sex': int(row['Sex']),
                'Population': row['Population'],
                'Superpopulation': row['Superpopulation']
            }
            
            # Obter likelihood FROG se disponível
            frog_likelihood = None
            frog_pop_names = None
            if frog_data and frog_data.has_sample(sample_id):
                try:
                    frog_info = frog_data.get_individual_data(sample_id, full_frog_vector=True)
                    frog_likelihood = frog_info['likelihood']
                    frog_pop_names = frog_info['population_names']
                except Exception as e:
                    print(f"[WARN] Erro ao obter likelihood para {sample_id}: {e}")
            
            # Criar builder de metadados individuais
            ind_builder = IndividualDatasetBuilder(
                base_dir=output_dir,
                sample_id=sample_id,
                sample_info=sample_info,
                frog_likelihood=frog_likelihood,
                frog_population_names=frog_pop_names
            )
            ind_builder.create_structure()
            
            # Executar build_window_and_predict
            success, target_name = run_build_window_predict(sample_id, config, output_dir)
            
            if success:
                # Descobrir quais janelas foram criadas
                windows_dir = output_dir / "individuals" / sample_id / "windows"
                if windows_dir.exists():
                    window_names = [d.name for d in windows_dir.iterdir() if d.is_dir()]
                    
                    # Adicionar metadados de cada janela
                    for window_name in window_names:
                        # Tentar determinar informações da janela
                        mode = params.get('mode', 'gene')
                        
                        # Para simplificar, usaremos valores padrão
                        # Em produção, poderíamos parsear os arquivos para obter estas informações
                        ind_builder.add_window(
                            target_name=window_name,
                            window_type=mode,
                            chromosome='unknown',  # Seria ideal extrair do FASTA ou metadados
                            start=0,
                            end=1000000,
                            outputs=params.get('outputs', '').split(',') if params.get('outputs') else [],
                            ontologies=params.get('ontology', '').split(',') if params.get('ontology') else []
                        )
                
                # Salvar metadados individuais
                ind_builder.save_metadata()
                
                checkpoint_data['completed_samples'].append(sample_id)
                print(f"[INFO] ✓ Amostra {sample_id} processada com sucesso")
            else:
                checkpoint_data['failed_samples'].append(sample_id)
                print(f"[ERROR] ✗ Falha ao processar amostra {sample_id}")
            
            # Save checkpoint
            save_checkpoint(checkpoint_file, checkpoint_data)
        
        print(f"\n[INFO] Processamento concluído!")
        print(f"  • Sucesso: {len(checkpoint_data['completed_samples'])}")
        print(f"  • Falha: {len(checkpoint_data['failed_samples'])}")
    
    # Step 5: Generate report
    if steps.get('generate_report', False):
        print("\n" + "="*80)
        print("PASSO 5: GERAÇÃO DE RELATÓRIO")
        print("="*80)
        
        checkpoint_file = output_dir / config['pipeline']['checkpoint_file']
        checkpoint_data = load_checkpoint(checkpoint_file)
        
        generate_summary_report(config, checkpoint_data, output_dir)
    
    # Step 6: Generate dataset metadata (PyTorch Dataset)
    if steps.get('generate_dataset_metadata', False):
        print("\n" + "="*80)
        print("PASSO 6: GERAÇÃO DE METADADOS DO DATASET PYTORCH")
        print("="*80)
        
        try:
            dataset_name = config['project'].get('name', 'non_longevous_dataset')
            
            metadata_builder = DatasetMetadataBuilder(
                dataset_dir=output_dir,
                dataset_name=dataset_name
            )
            
            # Escanear todos os indivíduos
            metadata_builder.scan_individuals()
            
            # Salvar metadados globais
            metadata_builder.save_metadata()
            
            # Imprimir sumário
            metadata_builder.print_summary()
            
            print(f"[INFO] ✓ Metadados do Dataset PyTorch gerados com sucesso!")
            print(f"[INFO] Arquivo: {output_dir / 'dataset_metadata.json'}")
            
        except Exception as e:
            print(f"[ERROR] Erro ao gerar metadados do dataset: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n[DONE] Pipeline concluído!")


if __name__ == "__main__":
    main()

