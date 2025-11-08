# === Configuração ===
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
VEPCACHE="/dados/vep_cache"
mkdir -p "$VEPCACHE"

# Checagem básica do Perl DBI (VEP precisa)
perl -MDBI -e 1 || { echo "Perl DBI ausente no env. Ative o 'genomics' certo ou reinstale ensembl-vep.";}

# Mostra de onde vem o vep (só informativo)
echo "vep em: $(command -v vep)"

# vep_install. Se perguntar "Do you wish to exit so you can get updates (y) or continue (n):" escolha n. Demora pra caramba... ===
vep_install -a cf -s "$SPECIES" -y "$ASSEMBLY" -c "$VEPCACHE" --NO_BIOPERL

# === Verificação rápida do cache baixado ===
echo "Conteúdo esperado:"
ls -l "$VEPCACHE/$SPECIES" || true
ls -l "$VEPCACHE/$SPECIES"/*"${ASSEMBLY}"* 2>/dev/null || true

