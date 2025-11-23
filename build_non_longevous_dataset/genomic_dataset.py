#!/usr/bin/env python3
"""
genomic_dataset.py

Dataset PyTorch para dados genômicos de longevidade.
Cada sample representa um indivíduo com múltiplas janelas genômicas.

Author: ChatGPT (for Alberto)
Created: 2025-11-11
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings


try:
    import torch
    from torch.utils.data import Dataset
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    warnings.warn(
        "PyTorch não está instalado. GenomicLongevityDataset não estará disponível. "
        "Instale com: pip install torch"
    )
    # Definir classe base mock para evitar erros
    class Dataset:
        pass


class GenomicLongevityDataset(Dataset):
    """
    Dataset PyTorch para dados genômicos de longevidade.
    
    Cada sample corresponde a um indivíduo com múltiplas janelas genômicas (1Mb),
    predições AlphaGenome, e metadados de ancestralidade.
    
    Attributes:
        dataset_dir: Diretório raiz do dataset
        load_predictions: Se True, carrega predições AlphaGenome
        load_sequences: Se True, carrega sequências FASTA
        transform: Transformação opcional a aplicar aos dados
        individuals: Lista de IDs de indivíduos
        dataset_metadata: Metadados globais do dataset
    
    Example:
        >>> from genomic_dataset import GenomicLongevityDataset
        >>> dataset = GenomicLongevityDataset(
        ...     dataset_dir='non_longevous_results',
        ...     load_predictions=True,
        ...     load_sequences=True
        ... )
        >>> print(f"Total: {len(dataset)} indivíduos")
        >>> input_data, output_data = dataset[0]
        >>> print(f"Sample: {output_data['sample_id']}")
    """
    
    def __init__(
        self,
        dataset_dir: Path,
        load_predictions: bool = True,
        load_sequences: bool = True,
        transform: Optional[Any] = None,
        cache_sequences: bool = False
    ):
        """
        Inicializa o dataset.
        
        Args:
            dataset_dir: Diretório contendo individuals/ e dataset_metadata.json
            load_predictions: Se True, carrega predições AlphaGenome (.npz)
            load_sequences: Se True, carrega sequências FASTA
            transform: Transformação opcional (callable) para aplicar aos dados
            cache_sequences: Se True, mantém sequências em memória (usar com cuidado)
        
        Raises:
            FileNotFoundError: Se dataset_dir não existir ou estrutura inválida
            ValueError: Se não houver indivíduos válidos
        """
        if not TORCH_AVAILABLE:
            raise ImportError(
                "PyTorch não está instalado. Instale com: pip install torch"
            )
        
        self.dataset_dir = Path(dataset_dir)
        self.load_predictions = load_predictions
        self.load_sequences = load_sequences
        self.transform = transform
        self.cache_sequences = cache_sequences
        
        # Validar estrutura
        self._validate_structure()
        
        # Carregar metadados globais
        self.dataset_metadata = self._load_dataset_metadata()
        
        # Listar indivíduos
        self.individuals = self.dataset_metadata.get('individuals', [])
        
        if len(self.individuals) == 0:
            raise ValueError(f"Nenhum indivíduo encontrado em {self.dataset_dir}")
        
        # Cache de sequências (se habilitado)
        self._sequence_cache = {} if cache_sequences else None
        
        print(f"[INFO] GenomicLongevityDataset inicializado:")
        print(f"  • Dataset: {self.dataset_metadata['dataset_name']}")
        print(f"  • Indivíduos: {len(self.individuals)}")
        print(f"  • Load predictions: {load_predictions}")
        print(f"  • Load sequences: {load_sequences}")
        print(f"  • Cache sequences: {cache_sequences}")
    
    def _validate_structure(self) -> None:
        """
        Valida estrutura básica do dataset.
        
        Raises:
            FileNotFoundError: Se estrutura for inválida
        """
        if not self.dataset_dir.exists():
            raise FileNotFoundError(f"Dataset directory não encontrado: {self.dataset_dir}")
        
        individuals_dir = self.dataset_dir / "individuals"
        if not individuals_dir.exists():
            raise FileNotFoundError(
                f"Diretório 'individuals' não encontrado em {self.dataset_dir}"
            )
        
        metadata_file = self.dataset_dir / "dataset_metadata.json"
        if not metadata_file.exists():
            raise FileNotFoundError(
                f"Arquivo dataset_metadata.json não encontrado em {self.dataset_dir}"
            )
    
    def _load_dataset_metadata(self) -> Dict:
        """
        Carrega metadados globais do dataset.
        
        Returns:
            Dicionário com metadados
        """
        metadata_file = self.dataset_dir / "dataset_metadata.json"
        with open(metadata_file, 'r') as f:
            return json.load(f)
    
    def _load_individual_metadata(self, sample_id: str) -> Dict:
        """
        Carrega metadados de um indivíduo.
        
        Args:
            sample_id: ID do indivíduo
            
        Returns:
            Dicionário com metadados do indivíduo
        """
        metadata_file = self.dataset_dir / "individuals" / sample_id / "individual_metadata.json"
        
        if not metadata_file.exists():
            raise FileNotFoundError(
                f"Metadados não encontrados para {sample_id}: {metadata_file}"
            )
        
        with open(metadata_file, 'r') as f:
            return json.load(f)
    
    def _load_fasta_sequence(self, fasta_path: Path) -> str:
        """
        Carrega sequência de um arquivo FASTA.
        
        Args:
            fasta_path: Caminho para arquivo .fa
            
        Returns:
            String com a sequência (sem header)
        """
        # Verificar cache
        if self._sequence_cache is not None and str(fasta_path) in self._sequence_cache:
            return self._sequence_cache[str(fasta_path)]
        
        if not fasta_path.exists():
            raise FileNotFoundError(f"Arquivo FASTA não encontrado: {fasta_path}")
        
        sequence_lines = []
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    continue  # Pular header
                sequence_lines.append(line)
        
        sequence = ''.join(sequence_lines)
        
        # Adicionar ao cache se habilitado
        if self._sequence_cache is not None:
            self._sequence_cache[str(fasta_path)] = sequence
        
        return sequence
    
    def _load_predictions(self, predictions_dir: Path) -> Dict[str, np.ndarray]:
        """
        Carrega predições AlphaGenome de um diretório.
        
        Args:
            predictions_dir: Diretório contendo arquivos .npz
            
        Returns:
            Dicionário com {output_type: array}
            Ex: {'rna_seq': array(...), 'atac': array(...)}
            
            Arrays podem ser 1D (1 track) ou 2D (múltiplas tracks como colunas).
        """
        if not predictions_dir.exists():
            return {}
        
        predictions = {}
        
        # Listar arquivos .npz
        npz_files = list(predictions_dir.glob("*.npz"))
        
        for npz_file in npz_files:
            output_type = npz_file.stem  # Nome sem extensão
            
            try:
                data = np.load(npz_file)
                
                # Tentar carregar 'values' se existir, senão pegar primeiro array
                if 'values' in data:
                    predictions[output_type] = data['values']
                else:
                    # Pegar primeiro array disponível
                    keys = list(data.keys())
                    if keys:
                        predictions[output_type] = data[keys[0]]
                
            except Exception as e:
                warnings.warn(f"Erro ao carregar {npz_file}: {e}")
                continue
        
        return predictions
    
    def _load_window_data(
        self,
        sample_id: str,
        window_name: str
    ) -> Dict:
        """
        Carrega todos os dados de uma janela.
        
        Args:
            sample_id: ID do indivíduo
            window_name: Nome da janela (gene/SNP)
            
        Returns:
            Dicionário com dados da janela
        """
        window_dir = self.dataset_dir / "individuals" / sample_id / "windows" / window_name
        
        window_data = {}
        
        # Carregar sequências (se solicitado)
        if self.load_sequences:
            ref_fasta = window_dir / "ref.window.fa"
            h1_fasta = window_dir / f"{sample_id}.H1.window.fixed.fa"
            h2_fasta = window_dir / f"{sample_id}.H2.window.fixed.fa"
            
            try:
                window_data['ref_sequence'] = self._load_fasta_sequence(ref_fasta)
            except FileNotFoundError:
                window_data['ref_sequence'] = None
            
            try:
                window_data['h1_sequence'] = self._load_fasta_sequence(h1_fasta)
            except FileNotFoundError:
                window_data['h1_sequence'] = None
            
            try:
                window_data['h2_sequence'] = self._load_fasta_sequence(h2_fasta)
            except FileNotFoundError:
                window_data['h2_sequence'] = None
        
        # Carregar predições (se solicitado)
        if self.load_predictions:
            predictions_h1_dir = window_dir / "predictions_H1"
            predictions_h2_dir = window_dir / "predictions_H2"
            
            window_data['predictions_h1'] = self._load_predictions(predictions_h1_dir)
            window_data['predictions_h2'] = self._load_predictions(predictions_h2_dir)
        
        return window_data
    
    def __len__(self) -> int:
        """
        Retorna número de indivíduos no dataset.
        
        Returns:
            Número de samples
        """
        return len(self.individuals)
    
    def __getitem__(self, idx: int) -> Tuple[Dict, Dict]:
        """
        Retorna dados de entrada e saída para um indivíduo.
        
        Args:
            idx: Índice do indivíduo (0 a len-1)
            
        Returns:
            Tupla (input_data, output_data):
            
            input_data (dict): {
                'windows': {
                    'CYP2B6': {
                        'ref_sequence': str (ou None),
                        'h1_sequence': str (ou None),
                        'h2_sequence': str (ou None),
                        'predictions_h1': {
                            'rna_seq': np.ndarray,
                            'atac': np.ndarray,
                            ...
                        },
                        'predictions_h2': {...}
                    },
                    'rs10497191': {...},
                    ...
                }
            }
            
            output_data (dict): {
                'sample_id': str,
                'longevity': int (0 para não longevos),
                'sex': int (1=Male, 2=Female),
                'population': str,
                'superpopulation': str,
                'frog_likelihood': np.ndarray (shape ~150),
                'frog_population_names': List[str]
            }
        """
        sample_id = self.individuals[idx]
        
        # Carregar metadados do indivíduo
        individual_metadata = self._load_individual_metadata(sample_id)
        
        # Preparar dados de saída (labels)
        output_data = {
            'sample_id': sample_id,
            'longevity': 0 if not individual_metadata.get('longevity', False) else 1,
            'sex': individual_metadata.get('sex', 0),
            'population': individual_metadata.get('population', ''),
            'superpopulation': individual_metadata.get('superpopulation', ''),
        }
        
        # Adicionar likelihood FROG se disponível
        if 'frog_likelihood' in individual_metadata:
            output_data['frog_likelihood'] = np.array(individual_metadata['frog_likelihood'])
            output_data['frog_population_names'] = individual_metadata.get('frog_population_names', [])
        else:
            output_data['frog_likelihood'] = None
            output_data['frog_population_names'] = []
        
        # Preparar dados de entrada
        input_data = {'windows': {}}
        
        # Carregar dados de cada janela
        windows = individual_metadata.get('windows', [])
        
        for window_name in windows:
            try:
                window_data = self._load_window_data(sample_id, window_name)
                input_data['windows'][window_name] = window_data
            except Exception as e:
                warnings.warn(f"Erro ao carregar janela {window_name} para {sample_id}: {e}")
                continue
        
        # Aplicar transformação se fornecida
        if self.transform is not None:
            input_data, output_data = self.transform(input_data, output_data)
        
        return input_data, output_data
    
    def get_sample_by_id(self, sample_id: str) -> Tuple[Dict, Dict]:
        """
        Retorna dados para um sample ID específico.
        
        Args:
            sample_id: ID do indivíduo
            
        Returns:
            Tupla (input_data, output_data)
            
        Raises:
            ValueError: Se sample_id não existir
        """
        if sample_id not in self.individuals:
            raise ValueError(f"Sample ID '{sample_id}' não encontrado no dataset")
        
        idx = self.individuals.index(sample_id)
        return self[idx]
    
    def get_summary(self) -> Dict:
        """
        Retorna sumário do dataset.
        
        Returns:
            Dicionário com estatísticas
        """
        return {
            'dataset_name': self.dataset_metadata['dataset_name'],
            'total_individuals': len(self.individuals),
            'populations': self.dataset_metadata.get('population_distribution', {}),
            'superpopulations': self.dataset_metadata.get('superpopulation_distribution', {}),
            'window_size': self.dataset_metadata.get('window_size', 1000000),
            'alphagenome_outputs': self.dataset_metadata.get('alphagenome_outputs', []),
            'frog_population_count': self.dataset_metadata.get('frog_population_count', 0)
        }
    
    def print_summary(self) -> None:
        """
        Imprime sumário formatado do dataset.
        """
        summary = self.get_summary()
        
        print("\n" + "="*80)
        print("DATASET PYTORCH - SUMÁRIO")
        print("="*80)
        print(f"Nome: {summary['dataset_name']}")
        print(f"Total de indivíduos: {summary['total_individuals']}")
        print(f"Tamanho da janela: {summary['window_size']:,} bp")
        
        print(f"\nSuperpopulações:")
        for superpop, count in sorted(summary['superpopulations'].items()):
            print(f"  • {superpop}: {count} indivíduos")
        
        print(f"\nAlphaGenome outputs disponíveis:")
        for output in summary['alphagenome_outputs']:
            print(f"  • {output}")
        
        print(f"\nPopulações FROG (likelihood): {summary['frog_population_count']}")
        print("="*80 + "\n")
    
    def clear_cache(self) -> None:
        """
        Limpa cache de sequências (se habilitado).
        """
        if self._sequence_cache is not None:
            self._sequence_cache.clear()
            print(f"[INFO] Cache de sequências limpo")


def collate_fn_variable_windows(batch: List[Tuple[Dict, Dict]]) -> Tuple[List[Dict], Dict]:
    """
    Collate function customizada para batch de indivíduos com número variável de janelas.
    
    Esta função NÃO empilha os dados em tensores (devido a tamanhos variáveis),
    apenas agrupa em listas. Use esta função com DataLoader.
    
    Args:
        batch: Lista de tuplas (input_data, output_data)
        
    Returns:
        Tupla (batch_input, batch_output):
        - batch_input: Lista de dicts de input
        - batch_output: Dict com listas de cada campo de output
    
    Example:
        >>> from torch.utils.data import DataLoader
        >>> dataloader = DataLoader(
        ...     dataset,
        ...     batch_size=4,
        ...     collate_fn=collate_fn_variable_windows
        ... )
    """
    batch_input = []
    batch_output = {
        'sample_id': [],
        'longevity': [],
        'sex': [],
        'population': [],
        'superpopulation': [],
        'frog_likelihood': []
    }
    
    for input_data, output_data in batch:
        batch_input.append(input_data)
        
        for key in batch_output.keys():
            batch_output[key].append(output_data.get(key))
    
    # Converter listas numéricas para tensors se possível
    if TORCH_AVAILABLE:
        batch_output['longevity'] = torch.tensor(batch_output['longevity'])
        batch_output['sex'] = torch.tensor(batch_output['sex'])
        
        # frog_likelihood: empilhar se todos tiverem mesmo tamanho
        likelihoods = [x for x in batch_output['frog_likelihood'] if x is not None]
        if likelihoods and all(len(x) == len(likelihoods[0]) for x in likelihoods):
            batch_output['frog_likelihood'] = torch.tensor(np.array(likelihoods))
    
    return batch_input, batch_output


# Exemplo de uso
if __name__ == "__main__":
    import sys
    
    print("="*80)
    print("TESTE: GenomicLongevityDataset")
    print("="*80)
    
    # Verificar se PyTorch está disponível
    if not TORCH_AVAILABLE:
        print("[ERROR] PyTorch não está instalado. Instale com: pip install torch")
        sys.exit(1)
    
    # Path de exemplo (ajuste conforme necessário)
    dataset_dir = Path("non_longevous_results")
    
    if not dataset_dir.exists():
        print(f"[ERROR] Dataset não encontrado: {dataset_dir}")
        print("[INFO] Execute build_non_longevous_dataset.py primeiro para criar o dataset.")
        sys.exit(1)
    
    try:
        # Carregar dataset
        dataset = GenomicLongevityDataset(
            dataset_dir=dataset_dir,
            load_predictions=True,
            load_sequences=False,  # Não carregar sequências para teste rápido
            cache_sequences=False
        )
        
        # Imprimir sumário
        dataset.print_summary()
        
        # Testar acesso a um sample
        if len(dataset) > 0:
            print(f"[TEST] Acessando primeiro sample...")
            input_data, output_data = dataset[0]
            
            print(f"\n  Sample ID: {output_data['sample_id']}")
            print(f"  População: {output_data['population']}")
            print(f"  Superpopulação: {output_data['superpopulation']}")
            print(f"  Janelas disponíveis: {list(input_data['windows'].keys())}")
            
            # Verificar predições
            for window_name, window_data in input_data['windows'].items():
                if 'predictions_h1' in window_data:
                    pred_h1 = window_data['predictions_h1']
                    print(f"\n  Janela '{window_name}' - Predições H1:")
                    for output_type, array in pred_h1.items():
                        print(f"    • {output_type}: shape={array.shape}, dtype={array.dtype}")
        
        print("\n[DONE] Dataset funcionando corretamente!")
        
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

