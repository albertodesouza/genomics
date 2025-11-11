#!/usr/bin/env python3
"""
frog_ancestry_parser.py

Parser para arquivos de likelihood do FROGAncestryCalc.
Carrega e processa dados de ancestralidade para uso no Dataset PyTorch.

Author: ChatGPT (for Alberto)
Created: 2025-11-11
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional, List


def load_frog_likelihoods(filepath: Path) -> pd.DataFrame:
    """
    Carrega o arquivo de likelihood do FROGAncestryCalc.
    
    Args:
        filepath: Caminho para whole_1000genomes_55aisnps_likelihood.txt
        
    Returns:
        DataFrame com colunas: Individual, SNP_Count, e uma coluna por população (~150)
        
    Example:
        >>> df = load_frog_likelihoods(Path("FROGAncestryCalc/output/whole_1000genomes_55aisnps_likelihood.txt"))
        >>> print(df.shape)
        (3202, 152)  # 3202 indivíduos, 152 colunas (Individual + SNP_Count + 150 populações)
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Arquivo de likelihood não encontrado: {filepath}")
    
    # Carregar arquivo tab-delimited
    df = pd.read_csv(filepath, sep='\t')
    
    print(f"[INFO] FROGAncestryCalc likelihoods carregado: {len(df)} indivíduos, {len(df.columns)-2} populações")
    
    return df


def get_population_columns(df: pd.DataFrame) -> List[str]:
    """
    Retorna lista de colunas de população (exclui Individual e SNP_Count).
    
    Args:
        df: DataFrame de likelihoods
        
    Returns:
        Lista com nomes das populações
    """
    # Primeiras duas colunas são Individual e SNP_Count
    return df.columns[2:].tolist()


def get_individual_likelihood(sample_id: str, df: pd.DataFrame) -> Tuple[np.ndarray, List[str]]:
    """
    Retorna vetor de likelihood para um indivíduo específico.
    
    Args:
        sample_id: ID do indivíduo (ex: "HG01879")
        df: DataFrame de likelihoods
        
    Returns:
        Tupla com:
        - np.ndarray com likelihood values (~150 valores)
        - Lista com nomes das populações correspondentes
        
    Raises:
        ValueError: Se sample_id não for encontrado
        
    Example:
        >>> likelihood, pop_names = get_individual_likelihood("HG01879", df)
        >>> print(likelihood.shape)
        (150,)
        >>> print(pop_names[:3])
        ['Taiwanese Han', 'Ami', 'Hakka']
    """
    if sample_id not in df['Individual'].values:
        raise ValueError(f"Sample ID '{sample_id}' não encontrado no arquivo de likelihoods")
    
    # Filtrar linha do indivíduo
    row = df[df['Individual'] == sample_id].iloc[0]
    
    # Extrair colunas de população (pular Individual e SNP_Count)
    population_cols = get_population_columns(df)
    likelihood_values = row[population_cols].values.astype(float)
    
    return likelihood_values, population_cols


def load_population_mapping(filepath: Path) -> Dict[str, Dict[str, str]]:
    """
    Carrega mapeamento de populações FROGAncestryCalc → 1000 Genomes.
    
    Args:
        filepath: Caminho para population_mapping_1000genomes.csv
        
    Returns:
        Dicionário mapeando população FROG para informações 1000G:
        {
            'Luhya(LWK)': {
                'code_1000g': 'LWK',
                'full_name_1000g': 'Luhya in Webuye Kenya',
                'superpopulation_code': 'AFR',
                'superpopulation_name': 'African'
            },
            ...
        }
        
    Example:
        >>> mapping = load_population_mapping(Path("FROGAncestryCalc/population_mapping_1000genomes.csv"))
        >>> print(mapping['Yoruba(YRI)']['code_1000g'])
        'YRI'
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Arquivo de mapeamento não encontrado: {filepath}")
    
    # Carregar CSV
    df = pd.read_csv(filepath)
    
    # Criar dicionário de mapeamento
    mapping = {}
    for _, row in df.iterrows():
        frog_pop = row['FROGAncestryCalc_Population']
        mapping[frog_pop] = {
            'code_1000g': row['Pop_Code_1000G'],
            'full_name_1000g': row['Full_Name_1000G'],
            'superpopulation_code': row['Superpopulation_1000G'],
            'superpopulation_name': row['Superpopulation_Full_Name']
        }
    
    print(f"[INFO] Mapeamento de populações carregado: {len(mapping)} populações mapeadas")
    
    return mapping


def get_1000g_populations_from_frog(
    frog_population_names: List[str],
    mapping: Dict[str, Dict[str, str]]
) -> List[Optional[str]]:
    """
    Converte nomes de população FROG para códigos 1000 Genomes.
    
    Args:
        frog_population_names: Lista de nomes de população do FROG
        mapping: Dicionário de mapeamento
        
    Returns:
        Lista com códigos 1000G correspondentes (None para populações não mapeadas)
        
    Example:
        >>> frog_pops = ['Yoruba(YRI)', 'British(GBR)', 'Unknown Pop']
        >>> codes = get_1000g_populations_from_frog(frog_pops, mapping)
        >>> print(codes)
        ['YRI', 'GBR', None]
    """
    codes = []
    for frog_pop in frog_population_names:
        if frog_pop in mapping:
            codes.append(mapping[frog_pop]['code_1000g'])
        else:
            codes.append(None)
    
    return codes


def filter_likelihood_to_1000g_populations(
    likelihood: np.ndarray,
    frog_population_names: List[str],
    mapping: Dict[str, Dict[str, str]]
) -> Tuple[np.ndarray, List[str]]:
    """
    Filtra vetor de likelihood para incluir apenas populações mapeadas no 1000 Genomes.
    
    Args:
        likelihood: Vetor de likelihood completo (~150 valores)
        frog_population_names: Lista de nomes de população correspondentes
        mapping: Dicionário de mapeamento
        
    Returns:
        Tupla com:
        - np.ndarray filtrado (apenas populações 1000G)
        - Lista com códigos 1000G correspondentes
        
    Example:
        >>> likelihood_filtered, codes_1000g = filter_likelihood_to_1000g_populations(
        ...     likelihood, pop_names, mapping
        ... )
        >>> print(len(codes_1000g))
        26  # Apenas as 26 populações do 1000 Genomes
    """
    if len(likelihood) != len(frog_population_names):
        raise ValueError(
            f"Tamanho de likelihood ({len(likelihood)}) não corresponde "
            f"ao número de populações ({len(frog_population_names)})"
        )
    
    filtered_likelihood = []
    filtered_codes = []
    
    for i, frog_pop in enumerate(frog_population_names):
        if frog_pop in mapping:
            filtered_likelihood.append(likelihood[i])
            filtered_codes.append(mapping[frog_pop]['code_1000g'])
    
    return np.array(filtered_likelihood), filtered_codes


class FROGAncestryData:
    """
    Classe helper para gerenciar dados de ancestralidade do FROGAncestryCalc.
    
    Attributes:
        likelihood_df: DataFrame com likelihoods
        population_mapping: Dicionário de mapeamento FROG → 1000G
        population_names: Lista de nomes de populações FROG
    """
    
    def __init__(
        self,
        likelihood_file: Path,
        mapping_file: Path
    ):
        """
        Inicializa carregando arquivos de likelihood e mapeamento.
        
        Args:
            likelihood_file: Caminho para whole_1000genomes_55aisnps_likelihood.txt
            mapping_file: Caminho para population_mapping_1000genomes.csv
        """
        self.likelihood_df = load_frog_likelihoods(likelihood_file)
        self.population_mapping = load_population_mapping(mapping_file)
        self.population_names = get_population_columns(self.likelihood_df)
        
        print(f"[INFO] FROGAncestryData inicializado:")
        print(f"  • {len(self.likelihood_df)} indivíduos")
        print(f"  • {len(self.population_names)} populações FROG")
        print(f"  • {len(self.population_mapping)} populações mapeadas para 1000G")
    
    def get_individual_data(
        self,
        sample_id: str,
        full_frog_vector: bool = True
    ) -> Dict:
        """
        Retorna todos os dados de ancestralidade para um indivíduo.
        
        Args:
            sample_id: ID do indivíduo
            full_frog_vector: Se True, retorna vetor completo (~150); 
                             Se False, retorna apenas populações 1000G (~26)
        
        Returns:
            Dicionário com:
            {
                'sample_id': str,
                'likelihood': np.ndarray,
                'population_names': List[str],
                'snp_count': int
            }
        """
        likelihood, pop_names = get_individual_likelihood(sample_id, self.likelihood_df)
        
        # Obter SNP count
        row = self.likelihood_df[self.likelihood_df['Individual'] == sample_id].iloc[0]
        snp_count = int(row['SNP_Count'])
        
        if not full_frog_vector:
            # Filtrar apenas populações 1000G
            likelihood, pop_names = filter_likelihood_to_1000g_populations(
                likelihood, pop_names, self.population_mapping
            )
        
        return {
            'sample_id': sample_id,
            'likelihood': likelihood,
            'population_names': pop_names,
            'snp_count': snp_count
        }
    
    def get_available_samples(self) -> List[str]:
        """
        Retorna lista de todos os sample IDs disponíveis.
        
        Returns:
            Lista de sample IDs
        """
        return self.likelihood_df['Individual'].tolist()
    
    def has_sample(self, sample_id: str) -> bool:
        """
        Verifica se um sample ID existe nos dados.
        
        Args:
            sample_id: ID do indivíduo
            
        Returns:
            True se existe, False caso contrário
        """
        return sample_id in self.likelihood_df['Individual'].values


# Exemplo de uso
if __name__ == "__main__":
    import sys
    
    # Caminhos relativos ao diretório do script
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    likelihood_file = project_root / "FROGAncestryCalc/output/whole_1000genomes_55aisnps_likelihood.txt"
    mapping_file = project_root / "FROGAncestryCalc/population_mapping_1000genomes.csv"
    
    if not likelihood_file.exists() or not mapping_file.exists():
        print("[ERROR] Arquivos do FROGAncestryCalc não encontrados.")
        print(f"  Likelihood: {likelihood_file}")
        print(f"  Mapping: {mapping_file}")
        sys.exit(1)
    
    # Teste básico
    print("="*80)
    print("TESTE: FROGAncestryCalc Parser")
    print("="*80)
    
    # Inicializar
    frog_data = FROGAncestryData(likelihood_file, mapping_file)
    
    # Testar com primeiro sample
    samples = frog_data.get_available_samples()
    test_sample = samples[0]
    
    print(f"\n[TEST] Sample: {test_sample}")
    
    # Vetor completo
    data_full = frog_data.get_individual_data(test_sample, full_frog_vector=True)
    print(f"  Vetor completo FROG: {len(data_full['likelihood'])} populações")
    print(f"  SNPs usados: {data_full['snp_count']}")
    print(f"  Top 3 likelihoods: {data_full['likelihood'][:3]}")
    
    # Apenas 1000G
    data_1000g = frog_data.get_individual_data(test_sample, full_frog_vector=False)
    print(f"  Vetor filtrado 1000G: {len(data_1000g['likelihood'])} populações")
    print(f"  Populações: {', '.join(data_1000g['population_names'][:5])}...")
    
    print("\n[DONE] Parser funcionando corretamente!")

