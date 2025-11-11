#!/usr/bin/env python3
"""
dataset_builder.py

Construtor de metadados para Dataset PyTorch de genômica.
Gerencia a estrutura de dados por indivíduo.

Author: ChatGPT (for Alberto)
Created: 2025-11-11
"""

import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import numpy as np


class IndividualDatasetBuilder:
    """
    Construtor de estrutura de dados para um indivíduo no dataset.
    
    Gerencia criação de diretórios, metadados, e organização de arquivos
    gerados por build_window_and_predict.py.
    """
    
    def __init__(
        self,
        base_dir: Path,
        sample_id: str,
        sample_info: Dict,
        frog_likelihood: Optional[np.ndarray] = None,
        frog_population_names: Optional[List[str]] = None
    ):
        """
        Inicializa builder para um indivíduo.
        
        Args:
            base_dir: Diretório base do dataset (ex: non_longevous_results)
            sample_id: ID do indivíduo (ex: "HG01879")
            sample_info: Dicionário com informações do CSV:
                {
                    'FamilyID': str,
                    'SampleID': str,
                    'Sex': int (1=Male, 2=Female),
                    'Population': str (ex: 'ACB'),
                    'Superpopulation': str (ex: 'AFR')
                }
            frog_likelihood: Vetor de likelihood do FROGAncestryCalc
            frog_population_names: Lista de nomes de populações correspondentes
        """
        self.base_dir = Path(base_dir)
        self.sample_id = sample_id
        self.sample_info = sample_info
        self.frog_likelihood = frog_likelihood
        self.frog_population_names = frog_population_names
        
        # Diretórios
        self.individuals_dir = self.base_dir / "individuals"
        self.individual_dir = self.individuals_dir / sample_id
        self.windows_dir = self.individual_dir / "windows"
        self.metadata_file = self.individual_dir / "individual_metadata.json"
        
        # Metadados
        self.metadata = {
            'sample_id': sample_id,
            'family_id': sample_info.get('FamilyID', '0'),
            'sex': int(sample_info.get('Sex', 0)),
            'population': sample_info.get('Population', ''),
            'superpopulation': sample_info.get('Superpopulation', ''),
            'longevity': False,  # Sempre False para non_longevous
            'windows': [],
            'window_metadata': {},
            'created_at': datetime.now().isoformat(),
            'last_updated': datetime.now().isoformat()
        }
        
        # Adicionar likelihood se fornecido
        if frog_likelihood is not None:
            self.metadata['frog_likelihood'] = frog_likelihood.tolist()
            if frog_population_names is not None:
                self.metadata['frog_population_names'] = frog_population_names
    
    def create_structure(self) -> None:
        """
        Cria estrutura de diretórios para o indivíduo.
        """
        self.individuals_dir.mkdir(parents=True, exist_ok=True)
        self.individual_dir.mkdir(exist_ok=True)
        self.windows_dir.mkdir(exist_ok=True)
        
        print(f"[INFO] Estrutura criada para {self.sample_id}: {self.individual_dir}")
    
    def get_window_dir(self, target_name: str) -> Path:
        """
        Retorna diretório para uma janela específica.
        
        Args:
            target_name: Nome do gene/SNP (ex: "CYP2B6", "rs10497191")
            
        Returns:
            Path para individuals/SAMPLE_ID/windows/TARGET_NAME/
        """
        return self.windows_dir / target_name
    
    def add_window(
        self,
        target_name: str,
        window_type: str,
        chromosome: str,
        start: int,
        end: int,
        outputs: List[str],
        ontologies: Optional[List[str]] = None
    ) -> None:
        """
        Adiciona informação sobre uma janela processada.
        
        Args:
            target_name: Nome do gene/SNP
            window_type: "gene" ou "snp"
            chromosome: Cromossomo (ex: "chr19")
            start: Posição inicial (0-based)
            end: Posição final (0-based)
            outputs: Lista de output types (ex: ["RNA_SEQ", "ATAC"])
            ontologies: Lista de ontologias usadas (ex: ["UBERON:0002107"])
        """
        if target_name not in self.metadata['windows']:
            self.metadata['windows'].append(target_name)
        
        self.metadata['window_metadata'][target_name] = {
            'type': window_type,
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'window_size': end - start,
            'outputs': outputs,
            'ontologies': ontologies or [],
            'added_at': datetime.now().isoformat()
        }
        
        self.metadata['last_updated'] = datetime.now().isoformat()
    
    def save_metadata(self) -> None:
        """
        Salva metadados do indivíduo em JSON.
        """
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
        
        print(f"[INFO] Metadados salvos: {self.metadata_file}")
    
    def load_metadata(self) -> Dict:
        """
        Carrega metadados existentes.
        
        Returns:
            Dicionário com metadados
            
        Raises:
            FileNotFoundError: Se arquivo não existir
        """
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Arquivo de metadados não encontrado: {self.metadata_file}")
        
        with open(self.metadata_file, 'r') as f:
            self.metadata = json.load(f)
        
        return self.metadata
    
    def metadata_exists(self) -> bool:
        """
        Verifica se arquivo de metadados já existe.
        
        Returns:
            True se existe, False caso contrário
        """
        return self.metadata_file.exists()
    
    def get_metadata(self) -> Dict:
        """
        Retorna metadados atuais.
        
        Returns:
            Dicionário com metadados
        """
        return self.metadata
    
    def validate_window_files(self, target_name: str) -> Dict[str, bool]:
        """
        Valida existência de arquivos esperados para uma janela.
        
        Args:
            target_name: Nome do gene/SNP
            
        Returns:
            Dicionário com status de cada tipo de arquivo:
            {
                'window_dir_exists': bool,
                'ref_fasta': bool,
                'h1_fasta': bool,
                'h2_fasta': bool,
                'predictions_h1_dir': bool,
                'predictions_h2_dir': bool
            }
        """
        window_dir = self.get_window_dir(target_name)
        
        validation = {
            'window_dir_exists': window_dir.exists(),
            'ref_fasta': (window_dir / 'ref.window.fa').exists(),
            'h1_fasta': (window_dir / f'{self.sample_id}.H1.window.fixed.fa').exists(),
            'h2_fasta': (window_dir / f'{self.sample_id}.H2.window.fixed.fa').exists(),
            'predictions_h1_dir': (window_dir / 'predictions_H1').exists(),
            'predictions_h2_dir': (window_dir / 'predictions_H2').exists()
        }
        
        return validation
    
    def get_summary(self) -> Dict:
        """
        Retorna sumário do indivíduo.
        
        Returns:
            Dicionário com informações resumidas
        """
        return {
            'sample_id': self.sample_id,
            'population': self.metadata.get('population', ''),
            'superpopulation': self.metadata.get('superpopulation', ''),
            'n_windows': len(self.metadata.get('windows', [])),
            'windows': self.metadata.get('windows', []),
            'has_frog_data': 'frog_likelihood' in self.metadata
        }


class DatasetMetadataBuilder:
    """
    Construtor de metadados globais do dataset.
    
    Agrega informações de todos os indivíduos.
    """
    
    def __init__(
        self,
        dataset_dir: Path,
        dataset_name: str = "non_longevous_1000g"
    ):
        """
        Inicializa builder de metadados globais.
        
        Args:
            dataset_dir: Diretório base do dataset
            dataset_name: Nome do dataset
        """
        self.dataset_dir = Path(dataset_dir)
        self.dataset_name = dataset_name
        self.metadata_file = self.dataset_dir / "dataset_metadata.json"
        
        self.metadata = {
            'dataset_name': dataset_name,
            'creation_date': datetime.now().isoformat(),
            'last_updated': datetime.now().isoformat(),
            'total_individuals': 0,
            'window_size': 1000000,
            'individuals': [],
            'population_distribution': {},
            'superpopulation_distribution': {},
            'window_types': set(),
            'alphagenome_outputs': set(),
            'ontologies': set(),
            'frog_population_count': 0
        }
    
    def scan_individuals(self) -> None:
        """
        Escaneia diretório individuals/ e coleta informações.
        """
        individuals_dir = self.dataset_dir / "individuals"
        
        if not individuals_dir.exists():
            print(f"[WARN] Diretório individuals não encontrado: {individuals_dir}")
            return
        
        # Listar diretórios de indivíduos
        individual_dirs = [d for d in individuals_dir.iterdir() if d.is_dir()]
        
        print(f"[INFO] Escaneando {len(individual_dirs)} indivíduos...")
        
        population_counts = {}
        superpopulation_counts = {}
        all_outputs = set()
        all_ontologies = set()
        all_window_types = set()
        frog_pop_count = 0
        
        valid_individuals = []
        
        for ind_dir in individual_dirs:
            metadata_file = ind_dir / "individual_metadata.json"
            
            if not metadata_file.exists():
                print(f"[WARN] Metadados não encontrados para {ind_dir.name}, pulando...")
                continue
            
            try:
                with open(metadata_file, 'r') as f:
                    ind_metadata = json.load(f)
                
                sample_id = ind_metadata['sample_id']
                valid_individuals.append(sample_id)
                
                # Contar populações
                pop = ind_metadata.get('population', 'Unknown')
                superpop = ind_metadata.get('superpopulation', 'Unknown')
                
                population_counts[pop] = population_counts.get(pop, 0) + 1
                superpopulation_counts[superpop] = superpopulation_counts.get(superpop, 0) + 1
                
                # Coletar tipos de janela e outputs
                for window_name, window_info in ind_metadata.get('window_metadata', {}).items():
                    all_window_types.add(window_info.get('type', 'unknown'))
                    all_outputs.update(window_info.get('outputs', []))
                    all_ontologies.update(window_info.get('ontologies', []))
                
                # Contar populações FROG
                if 'frog_population_names' in ind_metadata:
                    frog_pop_count = len(ind_metadata['frog_population_names'])
                
            except Exception as e:
                print(f"[ERROR] Erro ao processar {ind_dir.name}: {e}")
                continue
        
        # Atualizar metadados globais
        self.metadata['total_individuals'] = len(valid_individuals)
        self.metadata['individuals'] = sorted(valid_individuals)
        self.metadata['population_distribution'] = population_counts
        self.metadata['superpopulation_distribution'] = superpopulation_counts
        self.metadata['window_types'] = list(all_window_types)
        self.metadata['alphagenome_outputs'] = sorted(list(all_outputs))
        self.metadata['ontologies'] = sorted(list(all_ontologies))
        self.metadata['frog_population_count'] = frog_pop_count
        self.metadata['last_updated'] = datetime.now().isoformat()
        
        print(f"[INFO] Escaneamento completo:")
        print(f"  • Indivíduos válidos: {len(valid_individuals)}")
        print(f"  • Populações: {len(population_counts)}")
        print(f"  • Superpopulações: {len(superpopulation_counts)}")
        print(f"  • Tipos de janela: {all_window_types}")
        print(f"  • AlphaGenome outputs: {len(all_outputs)}")
    
    def save_metadata(self) -> None:
        """
        Salva metadados globais em JSON.
        """
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
        
        print(f"[INFO] Metadados globais salvos: {self.metadata_file}")
    
    def load_metadata(self) -> Dict:
        """
        Carrega metadados globais existentes.
        
        Returns:
            Dicionário com metadados
        """
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Arquivo de metadados não encontrado: {self.metadata_file}")
        
        with open(self.metadata_file, 'r') as f:
            self.metadata = json.load(f)
        
        return self.metadata
    
    def print_summary(self) -> None:
        """
        Imprime sumário do dataset.
        """
        print("\n" + "="*80)
        print("SUMÁRIO DO DATASET")
        print("="*80)
        print(f"Nome: {self.metadata['dataset_name']}")
        print(f"Data de criação: {self.metadata['creation_date']}")
        print(f"Total de indivíduos: {self.metadata['total_individuals']}")
        print(f"Tamanho da janela: {self.metadata['window_size']:,} bp")
        
        print(f"\nDistribuição por superpopulação:")
        for superpop, count in sorted(self.metadata['superpopulation_distribution'].items()):
            print(f"  • {superpop}: {count} indivíduos")
        
        print(f"\nDistribuição por população:")
        for pop, count in sorted(self.metadata['population_distribution'].items()):
            print(f"  • {pop}: {count} indivíduos")
        
        print(f"\nAlphaGenome outputs: {', '.join(self.metadata['alphagenome_outputs'])}")
        print(f"Ontologias usadas: {', '.join(self.metadata['ontologies'])}")
        print(f"Populações FROG: {self.metadata['frog_population_count']}")
        print("="*80 + "\n")


# Exemplo de uso
if __name__ == "__main__":
    import sys
    
    # Teste básico
    print("="*80)
    print("TESTE: Dataset Builder")
    print("="*80)
    
    # Criar builder para indivíduo de teste
    base_dir = Path("test_dataset")
    sample_info = {
        'FamilyID': 'BB01',
        'SampleID': 'TEST001',
        'Sex': 1,
        'Population': 'GBR',
        'Superpopulation': 'EUR'
    }
    
    frog_likelihood = np.random.random(150)
    frog_population_names = [f"Pop_{i}" for i in range(150)]
    
    builder = IndividualDatasetBuilder(
        base_dir=base_dir,
        sample_id='TEST001',
        sample_info=sample_info,
        frog_likelihood=frog_likelihood,
        frog_population_names=frog_population_names
    )
    
    # Criar estrutura
    builder.create_structure()
    
    # Adicionar janela
    builder.add_window(
        target_name='CYP2B6',
        window_type='gene',
        chromosome='chr19',
        start=41497445,
        end=42497445,
        outputs=['RNA_SEQ', 'ATAC'],
        ontologies=['UBERON:0002107']
    )
    
    # Salvar metadados
    builder.save_metadata()
    
    # Sumário
    summary = builder.get_summary()
    print(f"\n[TEST] Sumário do indivíduo:")
    print(f"  Sample: {summary['sample_id']}")
    print(f"  População: {summary['population']}")
    print(f"  Janelas: {summary['n_windows']}")
    print(f"  FROG data: {summary['has_frog_data']}")
    
    # Testar metadados globais
    dataset_builder = DatasetMetadataBuilder(base_dir)
    dataset_builder.scan_individuals()
    dataset_builder.save_metadata()
    dataset_builder.print_summary()
    
    # Cleanup
    if base_dir.exists():
        shutil.rmtree(base_dir)
        print(f"\n[CLEANUP] Diretório de teste removido: {base_dir}")
    
    print("\n[DONE] Builder funcionando corretamente!")

