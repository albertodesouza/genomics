from __future__ import annotations

from typing import Dict, Tuple

import torch
import torch.nn as nn
from rich.console import Console


console = Console()

class NNAncestryPredictor(nn.Module):
    """
    Rede neural totalmente conectada (NN) para predição de ancestralidade.
    
    Arquitetura:
    - Flatten (converte entrada 2D para 1D)
    - Camada de entrada (tamanho variável)
    - Camadas ocultas (configurável)
    - Camada de saída (softmax para classificação ou linear para regressão)
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(NNAncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)

        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        if len(GENE_ORDER) == 0 or input_shape[0] % len(GENE_ORDER) != 0:
            raise ValueError(
                f"Input shape incompatível com gene_order: {input_shape[0]} linhas para {len(GENE_ORDER)} genes"
            )
        tracks_per_gene = input_shape[0] // len(GENE_ORDER)
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair com base no número real de tracks por gene.
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        self.input_size = num_rows * effective_size  # Tamanho após flatten
        
        # Parâmetros de arquitetura
        hidden_layers = config['model']['hidden_layers']
        activation_type = config['model']['activation']
        dropout_rate = config['model']['dropout_rate']
        
        # Escolher função de ativação (para camadas intermediárias)
        if activation_type == 'relu':
            self.activation = nn.ReLU()
        elif activation_type == 'tanh':
            self.activation = nn.Tanh()
        elif activation_type == 'sigmoid':
            self.activation = nn.Sigmoid()
        else:
            raise ValueError(f"Ativação não suportada: {activation_type}")
        
        # Construir camadas
        # Arquitetura: camadas hidden (TODAS com ativação) + camada de saída (LINEAR antes do softmax)
        layers = []
        prev_size = self.input_size
        
        # TODAS as camadas hidden (com ativação configurada)
        for hidden_size in hidden_layers:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(self.activation)
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Camada de saída (sempre LINEAR, sem ativação - logits)
        layers.append(nn.Linear(prev_size, num_classes))
        
        self.network = nn.Sequential(*layers)
        
        # Inicializar pesos apropriadamente
        self._initialize_weights(activation_type)
        
        console.print(f"[green]Modelo NN criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {self.input_shape[0]} x {self.input_shape[1]}")
        console.print(f"  • Input size (após flatten): {self.input_size}")
        console.print(f"  • Hidden layers: {hidden_layers}")
        console.print(f"  • Ativação: {activation_type} (em todas as camadas hidden)")
        console.print(f"  • Arquitetura: flatten → camadas hidden → saída (logits)")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Flatten: [batch, num_rows_selected, effective_size] -> [batch, num_rows_selected * effective_size]
        x = x.view(x.size(0), -1)
        
        logits = self.network(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self, activation_type: str):
        """
        Inicializa pesos apropriadamente de acordo com tipo de ativação.
        
        Estratégia:
        - Camadas hidden (com ativação):
            * ReLU: He/Kaiming initialization (fan_in)
            * Tanh/Sigmoid: Xavier/Glorot initialization
        - Última camada (saída, linear antes do softmax): Xavier initialization
        - Bias: zeros (padrão recomendado)
        
        Args:
            activation_type: Tipo de função de ativação das camadas hidden
        """
        layer_count = 0
        total_layers = sum(1 for m in self.modules() if isinstance(m, nn.Linear))
        
        for m in self.modules():
            if isinstance(m, nn.Linear):
                layer_count += 1
                
                # Última camada (saída): sempre Xavier
                # (linear antes do softmax, sem ativação adicional)
                if layer_count == total_layers:
                    nn.init.xavier_normal_(m.weight)
                # Camadas hidden: depende da ativação configurada
                elif activation_type == 'relu':
                    # He initialization (Kaiming) para ReLU
                    # Usa fan_in para manter variância durante forward pass
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    # Xavier/Glorot initialization para tanh/sigmoid
                    nn.init.xavier_normal_(m.weight)
                else:
                    # Fallback: Xavier para qualquer outra ativação
                    nn.init.xavier_normal_(m.weight)
                
                # Bias: inicializar com zeros (padrão recomendado)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos inicializados:[/green]")
        console.print(f"  • Camadas hidden: {activation_type} initialization")
        console.print(f"  • Camada de saída: Xavier initialization")
        console.print(f"  • Bias: zeros")


class CNNAncestryPredictor(nn.Module):
    """
    Rede neural convolucional (CNN) para predição de ancestralidade.
    
    Arquitetura:
    - Input: [batch, num_rows, effective_size] como imagem com 1 canal
    - Conv2D: camada convolucional com kernel configurável
    - Ativação (ReLU/Tanh/Sigmoid)
    - MaxPool2D (opcional)
    - Flatten
    - Camadas fully connected (hidden_layers)
    - Output linear → Softmax (classificação)
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo CNN.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(CNNAncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)

        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        if len(GENE_ORDER) == 0 or input_shape[0] % len(GENE_ORDER) != 0:
            raise ValueError(
                f"Input shape incompatível com gene_order: {input_shape[0]} linhas para {len(GENE_ORDER)} genes"
            )
        tracks_per_gene = input_shape[0] // len(GENE_ORDER)
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair com base no número real de tracks por gene.
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        
        # Parâmetros CNN
        cnn_config = config['model']['cnn']
        kernel_size = tuple(cnn_config['kernel_size'])  # [height, width]
        num_filters = cnn_config['num_filters']
        
        # Stride pode ser escalar ou lista [vertical, horizontal]
        stride_config = cnn_config['stride']
        if isinstance(stride_config, list):
            stride = tuple(stride_config)
        else:
            stride = (stride_config, stride_config)  # Usar mesmo valor para ambas dimensões
        
        # Padding pode ser escalar ou lista [vertical, horizontal]
        padding_config = cnn_config['padding']
        if isinstance(padding_config, list):
            padding = tuple(padding_config)
        else:
            padding = padding_config  # PyTorch aceita int ou tuple
        
        pool_size = cnn_config.get('pool_size')  # Pode ser None
        
        # Parâmetros gerais
        activation_type = config['model']['activation']
        dropout_rate = config['model']['dropout_rate']
        hidden_layers = config['model']['hidden_layers']
        
        # Escolher função de ativação
        if activation_type == 'relu':
            self.activation = nn.ReLU()
        elif activation_type == 'tanh':
            self.activation = nn.Tanh()
        elif activation_type == 'sigmoid':
            self.activation = nn.Sigmoid()
        else:
            raise ValueError(f"Ativação não suportada: {activation_type}")
        
        # Camada convolucional
        # Input: [batch, 1, num_rows, effective_size]
        # Output: [batch, num_filters, out_height, out_width]
        self.conv1 = nn.Conv2d(
            in_channels=1,
            out_channels=num_filters,
            kernel_size=kernel_size,
            stride=stride,
            padding=padding
        )
        
        # Calcular dimensões após convolução (usando input_shape ajustado)
        num_rows, effective_size = self.input_shape
        
        # Extrair padding como tupla
        if isinstance(padding, int):
            pad_h, pad_w = padding, padding
        else:
            pad_h, pad_w = padding
        
        # Extrair stride como tupla
        if isinstance(stride, int):
            stride_h, stride_w = stride, stride
        else:
            stride_h, stride_w = stride
        
        # Fórmula: out_size = (in_size + 2*padding - kernel_size) / stride + 1
        conv_out_h = (num_rows + 2*pad_h - kernel_size[0]) // stride_h + 1
        conv_out_w = (effective_size + 2*pad_w - kernel_size[1]) // stride_w + 1
        
        # Pooling (opcional)
        if pool_size is not None:
            pool_size = tuple(pool_size)
            self.pool = nn.MaxPool2d(kernel_size=pool_size)
            pool_out_h = conv_out_h // pool_size[0]
            pool_out_w = conv_out_w // pool_size[1]
        else:
            self.pool = None
            pool_out_h = conv_out_h
            pool_out_w = conv_out_w
        
        # Tamanho após flatten
        flattened_size = num_filters * pool_out_h * pool_out_w
        
        # Camadas fully connected
        fc_layers = []
        prev_size = flattened_size
        
        for hidden_size in hidden_layers:
            fc_layers.append(nn.Linear(prev_size, hidden_size))
            fc_layers.append(self.activation)
            if dropout_rate > 0:
                fc_layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Camada de saída (linear, sem ativação - logits)
        fc_layers.append(nn.Linear(prev_size, num_classes))
        
        self.fc_network = nn.Sequential(*fc_layers)
        
        # Inicializar pesos
        self._initialize_weights(activation_type)
        
        console.print(f"[green]Modelo CNN criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {self.input_shape[0]} x {self.input_shape[1]} (1 canal)")
        console.print(f"  • Conv2D: {num_filters} filters, kernel={kernel_size}, stride={stride}, padding={padding}")
        console.print(f"  • Após Conv2D: {num_filters} x {conv_out_h} x {conv_out_w}")
        if pool_size:
            console.print(f"  • MaxPool2D: kernel={pool_size}")
            console.print(f"  • Após Pool: {num_filters} x {pool_out_h} x {pool_out_w}")
        console.print(f"  • Flatten size: {flattened_size}")
        console.print(f"  • FC hidden layers: {hidden_layers}")
        console.print(f"  • Ativação: {activation_type}")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Adicionar dimensão de canal: [batch, num_rows_selected, effective_size] -> [batch, 1, num_rows_selected, effective_size]
        x = x.unsqueeze(1)
        
        # Convolução + Ativação
        x = self.conv1(x)
        x = self.activation(x)
        
        # Pooling (se habilitado)
        if self.pool is not None:
            x = self.pool(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Camadas fully connected
        logits = self.fc_network(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self, activation_type: str):
        """
        Inicializa pesos apropriadamente de acordo com tipo de ativação.
        
        Estratégia:
        - Conv2D: He/Kaiming para ReLU, Xavier para Tanh/Sigmoid
        - FC hidden layers: He/Kaiming para ReLU, Xavier para Tanh/Sigmoid
        - Última camada FC (saída): Xavier initialization
        - Bias: zeros
        
        Args:
            activation_type: Tipo de função de ativação
        """
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                # Inicialização convolucional
                if activation_type == 'relu':
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.xavier_normal_(m.weight)
                
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            
            elif isinstance(m, nn.Linear):
                # Contar quantas camadas lineares existem
                linear_layers = [module for module in self.fc_network.modules() if isinstance(module, nn.Linear)]
                total_layers = len(linear_layers)
                layer_count = sum(1 for module in self.fc_network.modules() 
                                 if isinstance(module, nn.Linear) and id(module) <= id(m))
                
                # Última camada (saída): Xavier
                if layer_count == total_layers:
                    nn.init.xavier_normal_(m.weight)
                # Camadas hidden
                elif activation_type == 'relu':
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.xavier_normal_(m.weight)
                
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos CNN inicializados:[/green]")
        console.print(f"  • Conv2D: {activation_type} initialization")
        console.print(f"  • FC hidden layers: {activation_type} initialization")
        console.print(f"  • Camada de saída: Xavier initialization")
        console.print(f"  • Bias: zeros")


class CNN2AncestryPredictor(nn.Module):
    """
    Rede neural convolucional avançada (CNN2) para predição de ancestralidade.
    
    Arquitetura multi-layer com global pooling, mais adequada para dados de
    genes com múltiplas tracks (ex: 11 genes × 6 ontologias = 66 tracks).
    
    Arquitetura:
    - Input: [batch, 1, num_rows, effective_size]
    - Conv2D Stage 1: Agrupa tracks de genes (ex: 6→1 por gene)
    - Conv2D Stage 2: Processa ao longo das bases
    - Conv2D Stage 3: Features mais abstratas
    - AdaptiveAvgPool2d: Global pooling → vetor fixo
    - FC: Classificação final
    
    Vantagens sobre CNN simples:
    - Menos parâmetros (~100-200k vs 5M+)
    - Hierarquia de features (padrões locais → globais)
    - Global pooling = invariante ao tamanho da entrada
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo CNN2.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(CNN2AncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)

        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        if len(GENE_ORDER) == 0 or input_shape[0] % len(GENE_ORDER) != 0:
            raise ValueError(
                f"Input shape incompatível com gene_order: {input_shape[0]} linhas para {len(GENE_ORDER)} genes"
            )
        tracks_per_gene = input_shape[0] // len(GENE_ORDER)
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair com base no número real de tracks por gene.
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        
        # Ler parâmetros da configuração CNN2
        cnn2_config = config['model'].get('cnn2', {})
        
        # Stage 1: Primeira convolução (agrupa tracks)
        num_filters_s1 = cnn2_config.get('num_filters_stage1', 16)
        kernel_s1 = tuple(cnn2_config.get('kernel_stage1', [6, 32]))
        stride_s1 = tuple(cnn2_config.get('stride_stage1', [6, 32]))
        
        # Stage 2 e 3: Convoluções subsequentes
        num_filters_s2 = cnn2_config.get('num_filters_stage2', 32)
        num_filters_s3 = cnn2_config.get('num_filters_stage3', 64)
        kernel_s23 = tuple(cnn2_config.get('kernel_stages23', [1, 5]))
        stride_s23 = tuple(cnn2_config.get('stride_stages23', [1, 2]))
        padding_s23 = tuple(cnn2_config.get('padding_stages23', [0, 2]))
        
        # FC e dropout
        fc_hidden_size = cnn2_config.get('fc_hidden_size', 128)
        dropout_rate = config['model'].get('dropout_rate', 0.3)
        
        # Calcular dimensões após cada convolução
        # Stage 1
        conv1_h = (num_rows - kernel_s1[0]) // stride_s1[0] + 1
        conv1_w = (effective_size - kernel_s1[1]) // stride_s1[1] + 1
        
        # Stage 2
        conv2_h = (conv1_h + 2*padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv2_w = (conv1_w + 2*padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1
        
        # Stage 3
        conv3_h = (conv2_h + 2*padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv3_w = (conv2_w + 2*padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1

        if min(conv1_h, conv1_w, conv2_h, conv2_w, conv3_h, conv3_w) <= 0:
            raise ValueError(
                "CNN2 produziu dimensão inválida. "
                f"input=({num_rows}, {effective_size}), "
                f"stage1=({conv1_h}, {conv1_w}), "
                f"stage2=({conv2_h}, {conv2_w}), "
                f"stage3=({conv3_h}, {conv3_w}), "
                f"kernel_s1={kernel_s1}, stride_s1={stride_s1}, "
                f"kernel_s23={kernel_s23}, stride_s23={stride_s23}, padding_s23={padding_s23}"
            )
        
        # Bloco de features convolucionais
        self.features = nn.Sequential(
            # Stage 1: Agrupa tracks de genes
            nn.Conv2d(1, num_filters_s1, kernel_size=kernel_s1, stride=stride_s1),
            nn.ReLU(),
            
            # Stage 2: Processa ao longo das bases
            nn.Conv2d(num_filters_s1, num_filters_s2, 
                     kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU(),
            
            # Stage 3: Features mais abstratas
            nn.Conv2d(num_filters_s2, num_filters_s3,
                     kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU()
        )
        
        # Global pooling determinístico: usa MaxPool2d ou AvgPool2d com kernel fixo
        # Faz pooling apenas na dimensão da largura, mantendo a altura (genes)
        # NOTA: AdaptiveAvgPool2d não é determinístico em CUDA no backward pass
        pool_type = cnn2_config.get('global_pool_type', 'max')
        if pool_type == 'max':
            self.global_pool = nn.MaxPool2d(kernel_size=(1, conv3_w))
        else:
            self.global_pool = nn.AvgPool2d(kernel_size=(1, conv3_w))
        self.pool_type = pool_type
        
        # Dimensão após flatten
        flattened_size = num_filters_s3 * conv3_h
        
        # Classificador
        self.classifier = nn.Sequential(
            nn.Linear(flattened_size, fc_hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(fc_hidden_size, num_classes)
        )
        
        # Inicializar pesos
        self._initialize_weights()
        
        # Imprimir informações do modelo
        console.print(f"[green]Modelo CNN2 criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {num_rows} x {effective_size} (1 canal)")
        console.print(f"  • Stage 1: {num_filters_s1} filters, kernel={kernel_s1}, stride={stride_s1}")
        console.print(f"    → Output: {num_filters_s1} x {conv1_h} x {conv1_w}")
        console.print(f"  • Stage 2: {num_filters_s2} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23}")
        console.print(f"    → Output: {num_filters_s2} x {conv2_h} x {conv2_w}")
        console.print(f"  • Stage 3: {num_filters_s3} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23}")
        console.print(f"    → Output: {num_filters_s3} x {conv3_h} x {conv3_w}")
        console.print(f"  • Global Pool ({pool_type}): kernel=(1, {conv3_w}) → {num_filters_s3} x {conv3_h} x 1")
        console.print(f"  • Flatten size: {flattened_size}")
        console.print(f"  • FC: {flattened_size} → {fc_hidden_size} → {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Adicionar dimensão de canal: [batch, num_rows_selected, effective_size] -> [batch, 1, num_rows_selected, effective_size]
        x = x.unsqueeze(1)
        
        # Features convolucionais
        x = self.features(x)
        
        # Global pooling
        x = self.global_pool(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Classificador
        logits = self.classifier(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self):
        """
        Inicializa pesos usando He/Kaiming para ReLU.
        
        Estratégia:
        - Conv2D: Kaiming (He) initialization para ReLU
        - FC layers: Kaiming para camadas intermediárias, Xavier para saída
        - Bias: zeros
        """
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                # Última camada (saída) usa Xavier, outras usam Kaiming
                if m.out_features == self.num_classes:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos CNN2 inicializados:[/green]")
        console.print(f"  • Conv2D: Kaiming (He) initialization")
        console.print(f"  • FC hidden: Kaiming initialization")
        console.print(f"  • FC output: Xavier initialization")
        console.print(f"  • Bias: zeros")


