#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Teste simples de visualização do AlphaGenome
"""

import sys

try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome
    from alphagenome.visualization import plot_components
    import matplotlib.pyplot as plt
    
    API_KEY = sys.argv[1] if len(sys.argv) > 1 else "test"
    
    print("Criando cliente...")
    client = dna_client.create(API_KEY)
    
    # Criar sequência de teste (2048 bp)
    seq = "ATCG" * 512
    interval = genome.Interval(chromosome='chr1', start=1000000, end=1002048)
    ontology_terms = ['UBERON:0000955']
    
    print("Fazendo predição...")
    outputs = client.predict_interval(
        interval=interval,
        ontology_terms=ontology_terms,
        requested_outputs=[dna_client.OutputType.RNA_SEQ]
    )
    
    print("✓ Predição concluída")
    print(f"Tipo de rna_seq: {type(outputs.rna_seq)}")
    
    # Tentar plotar
    print("\nTeste 1: Usando Tracks com tdata")
    try:
        fig, ax = plt.subplots(figsize=(15, 5))
        plot_components.plot(
            [plot_components.Tracks(tdata={'RNA-seq': outputs.rna_seq})],
            interval=outputs.rna_seq.interval
        )
        plt.title("Test Tracks")
        plt.savefig("test_tracks.png", dpi=150, bbox_inches='tight')
        plt.close()
        print("✓ Tracks funcionou! Arquivo: test_tracks.png")
    except Exception as e:
        print(f"✗ Tracks falhou: {e}")
    
    # Teste 2: Sem wrapper
    print("\nTeste 2: Passando TrackData diretamente")
    try:
        fig, ax = plt.subplots(figsize=(15, 5))
        plot_components.plot(
            [outputs.rna_seq],
            interval=outputs.rna_seq.interval
        )
        plt.title("Test Direct TrackData")
        plt.savefig("test_direct.png", dpi=150, bbox_inches='tight')
        plt.close()
        print("✓ Direto funcionou! Arquivo: test_direct.png")
    except Exception as e:
        print(f"✗ Direto falhou: {e}")
    
    print("\nTestes concluídos!")
    
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

