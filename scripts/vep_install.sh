#!/usr/bin/env bash
# vep_install.sh — Instala o Ensembl VEP do GitHub e baixa o cache.
# Uso:  source scripts/vep_install.sh
#       source scripts/vep_install.sh [BRANCH]   (ex: release/114)
#
# Variáveis de ambiente opcionais (valores default mostrados):
#   VEP_BRANCH        release/115
#   VEP_SPECIES       homo_sapiens
#   VEP_ASSEMBLY      GRCh38
#   VEP_CACHE_DIR     $HOME/vep_cache
#   VEP_DIR           $HOME/ensembl-vep

# ---------------------------------------------------------------------------
# Configuração
# ---------------------------------------------------------------------------
VEP_BRANCH="${1:-${VEP_BRANCH:-release/115}}"
VEP_SPECIES="${VEP_SPECIES:-homo_sapiens}"
VEP_ASSEMBLY="${VEP_ASSEMBLY:-GRCh38}"
VEP_CACHE_DIR="${VEP_CACHE_DIR:-${HOME}/vep_cache}"
VEP_DIR="${VEP_DIR:-${HOME}/ensembl-vep}"

echo ""
echo "🧬 Instalador VEP — Ensembl GitHub"
echo "==================================="
echo "   Branch:   $VEP_BRANCH"
echo "   Espécie:  $VEP_SPECIES"
echo "   Assembly: $VEP_ASSEMBLY"
echo "   Cache:    $VEP_CACHE_DIR"
echo "   VEP dir:  $VEP_DIR"
echo ""

# ---------------------------------------------------------------------------
# Pré-requisitos
# ---------------------------------------------------------------------------
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "❌ Nenhum ambiente conda ativo."
    echo "💡 Execute: source scripts/start_genomics_universal.sh"
    return 1 2>/dev/null || exit 1
fi

if ! command -v git >/dev/null 2>&1; then
    echo "❌ git não encontrado."
    return 1 2>/dev/null || exit 1
fi

mkdir -p "$VEP_CACHE_DIR"

# ---------------------------------------------------------------------------
# 1. Clonar / atualizar repositório
# ---------------------------------------------------------------------------
echo "📥 1/4  Obtendo VEP do GitHub..."

if [ -d "$VEP_DIR/.git" ]; then
    echo "   Repositório já existe — atualizando..."
    git -C "$VEP_DIR" fetch --all --quiet
    git -C "$VEP_DIR" checkout "$VEP_BRANCH" --quiet 2>/dev/null || \
        git -C "$VEP_DIR" checkout "origin/$VEP_BRANCH" --quiet
    git -C "$VEP_DIR" pull --quiet 2>/dev/null || true
else
    [ -d "$VEP_DIR" ] && rm -rf "$VEP_DIR"
    git clone --depth 1 -b "$VEP_BRANCH" \
        https://github.com/Ensembl/ensembl-vep.git "$VEP_DIR"
fi

echo "   ✅ $(git -C "$VEP_DIR" log --oneline -1)"

# ---------------------------------------------------------------------------
# 2. Dependências Perl (conda — já instaladas por install_genomics_env.sh)
# ---------------------------------------------------------------------------
echo ""
echo "🐪 2/4  Verificando dependências Perl..."

_vep_check_perl() { perl -M"$1" -e 1 2>/dev/null; }

MISSING_PERL=()
for mod in DBI Archive::Zip JSON Try::Tiny; do
    _vep_check_perl "$mod" || MISSING_PERL+=("$mod")
done

if [ ${#MISSING_PERL[@]} -gt 0 ]; then
    echo "   Instalando módulos faltantes: ${MISSING_PERL[*]}"
    conda install -y -c conda-forge -c bioconda \
        perl-dbi perl-archive-zip perl-json perl-try-tiny 2>/dev/null || true
fi

echo "   ✅ Dependências OK"

# ---------------------------------------------------------------------------
# 3. Instalar APIs do VEP (não-interativo)
# ---------------------------------------------------------------------------
echo ""
echo "⚡ 3/4  Instalando APIs do VEP..."

chmod +x "$VEP_DIR/INSTALL.pl"

(cd "$VEP_DIR" && perl INSTALL.pl \
    --AUTO a \
    --NO_BIOPERL \
    --NO_HTSLIB \
    --NO_TEST \
    --NO_UPDATE)

echo "   ✅ APIs instaladas"

# ---------------------------------------------------------------------------
# 4. Baixar cache
# ---------------------------------------------------------------------------
echo ""
echo "💾 4/4  Baixando cache ($VEP_SPECIES / $VEP_ASSEMBLY)..."
echo "   Isso pode demorar 15–30 min na primeira execução."

(cd "$VEP_DIR" && perl INSTALL.pl \
    --AUTO cf \
    --SPECIES "$VEP_SPECIES" \
    --ASSEMBLY "$VEP_ASSEMBLY" \
    --CACHEDIR "$VEP_CACHE_DIR" \
    --NO_BIOPERL \
    --NO_HTSLIB \
    --NO_TEST \
    --NO_UPDATE)

echo "   ✅ Cache baixado"

# ---------------------------------------------------------------------------
# Configurar PATH no ambiente conda
# ---------------------------------------------------------------------------
echo ""
echo "🔗 Configurando PATH..."

export PATH="$VEP_DIR:$PATH"
export PERL5LIB="$VEP_DIR${PERL5LIB:+:$PERL5LIB}"

ACTIVATE_DIR="${CONDA_PREFIX}/etc/conda/activate.d"
if [ -d "${CONDA_PREFIX:-}" ]; then
    mkdir -p "$ACTIVATE_DIR"
    cat > "$ACTIVATE_DIR/vep.sh" <<VEPEOF
#!/bin/bash
export PATH="$VEP_DIR:\$PATH"
export PERL5LIB="$VEP_DIR\${PERL5LIB:+:\$PERL5LIB}"
VEPEOF
    echo "   ✅ PATH persistido em $ACTIVATE_DIR/vep.sh"
fi

# ---------------------------------------------------------------------------
# Verificação final
# ---------------------------------------------------------------------------
echo ""
echo "🧪 Verificação final..."

if command -v vep >/dev/null 2>&1; then
    echo "   ✅ vep encontrado: $(command -v vep)"
    VEP_VER=$(vep 2>&1 | grep -oP 'ensembl-vep\s*:\s*\K\S+' || echo "$VEP_BRANCH")
    echo "   ✅ Versão: $VEP_VER"
else
    echo "   ❌ vep não encontrado no PATH"
fi

if [ -d "$VEP_CACHE_DIR/$VEP_SPECIES" ]; then
    echo "   ✅ Cache: $(ls -d "$VEP_CACHE_DIR/$VEP_SPECIES"/*"$VEP_ASSEMBLY"* 2>/dev/null | head -1 || echo 'presente')"
else
    echo "   ⚠️  Cache não encontrado em $VEP_CACHE_DIR/$VEP_SPECIES"
fi

echo ""
echo "🎉 Instalação VEP concluída!"
echo "💡 Para usar:"
echo "   vep --cache --offline --species $VEP_SPECIES --assembly $VEP_ASSEMBLY --dir_cache $VEP_CACHE_DIR"
echo ""
