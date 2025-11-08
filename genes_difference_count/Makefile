# Makefile for genes_difference_count
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -mtune=native -flto -DNDEBUG
CXXFLAGS += -Wall -Wextra -pedantic
CXXFLAGS += -funroll-loops -ffast-math -finline-functions
CXXFLAGS += -fopenmp
TARGET = genes_difference_count
SOURCE = genes_difference_count.cpp

# Compilação otimizada para máxima performance
$(TARGET): $(SOURCE)
	@echo "🔨 Compilando versão ultra-otimizada..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)
	@echo "✅ Compilação concluída! Executável: $(TARGET)"

# Versão de debug
debug: CXXFLAGS = -std=c++17 -g -O0 -Wall -Wextra -pedantic -DDEBUG -fopenmp
debug: $(SOURCE)
	@echo "🐛 Compilando versão debug..."
	$(CXX) $(CXXFLAGS) -o $(TARGET)_debug $(SOURCE)
	@echo "✅ Versão debug compilada! Executável: $(TARGET)_debug"

# Limpeza
clean:
	@echo "🧹 Limpando arquivos compilados..."
	rm -f $(TARGET) $(TARGET)_debug
	@echo "✅ Limpeza concluída!"

# Teste rápido
test: $(TARGET)
	@echo "🧪 Executando teste..."
	./$(TARGET)

# Informações de compilação
info:
	@echo "📋 Informações de compilação:"
	@echo "   Compilador: $(CXX)"
	@echo "   Flags de otimização: $(CXXFLAGS)"
	@echo "   Target: $(TARGET)"
	@echo "   Source: $(SOURCE)"

# Instalar dependências (se necessário)
deps:
	@echo "📦 Verificando dependências..."
	@which g++ > /dev/null || (echo "❌ g++ não encontrado! Instale com: sudo apt install build-essential" && exit 1)
	@echo "✅ Todas as dependências estão OK!"

.PHONY: clean debug test info deps

# Explicação das otimizações:
# -O3: Máxima otimização
# -march=native: Otimizar para a CPU atual
# -mtune=native: Ajustar para a CPU atual  
# -flto: Link Time Optimization
# -DNDEBUG: Desabilitar asserts
# -funroll-loops: Desenrolar loops
# -ffast-math: Matemática rápida
# -finline-functions: Inlining agressivo
# -fopenmp: Habilitar OpenMP para paralelização
