#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Teste de visualização manual (sem plot_components)
"""

import sys
import numpy as np

try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome
    import matplotlib.pyplot as plt
    
    API_KEY = sys.argv[1] if len(sys.argv) > 1 else "test"
    
    print("Criando cliente...")
    client = dna_client.create(API_KEY)
    
    # Criar sequência de teste (2048 bp)
    interval = genome.Interval(chromosome='chr1', start=1000000, end=1002048)
    ontology_terms = ['UBERON:0000955']
    
    print("Fazendo predição...")
    outputs = client.predict_interval(
        interval=interval,
        ontology_terms=ontology_terms,
        requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.ATAC]
    )
    
    print("✓ Predição concluída\n")
    
    # Explorar TrackData
    rna_data = outputs.rna_seq
    print(f"Tipo: {type(rna_data)}")
    print(f"Atributos: {[a for a in dir(rna_data) if not a.startswith('_')]}")
    
    # Tentar acessar os dados
    print("\nTentando acessar dados...")
    try:
        # Muitos objetos TrackData tem .values ou .data
        if hasattr(rna_data, 'values'):
            print(f"✓ Tem .values: {type(rna_data.values)}")
            if hasattr(rna_data.values, 'shape'):
                print(f"  Shape: {rna_data.values.shape}")
        
        if hasattr(rna_data, 'data'):
            print(f"✓ Tem .data: {type(rna_data.data)}")
            if hasattr(rna_data.data, 'shape'):
                print(f"  Shape: {rna_data.data.shape}")
        
        # Tentar converter para numpy
        if hasattr(rna_data, '__array__'):
            print(f"✓ Pode ser convertido para numpy array")
            arr = np.array(rna_data)
            print(f"  Shape: {arr.shape}")
            
            # Tentar plotar
            fig, ax = plt.subplots(figsize=(15, 5))
            ax.plot(arr, linewidth=0.5)
            ax.set_title("RNA-seq - Plot Manual")
            ax.set_xlabel("Posição (bp)")
            ax.set_ylabel("Sinal")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig("manual_plot.png", dpi=150)
            plt.close()
            print("✓ Plot salvo em: manual_plot.png")
            
    except Exception as e:
        print(f"✗ Erro: {e}")
        import traceback
        traceback.print_exc()
    
except Exception as e:
    print(f"Erro fatal: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

