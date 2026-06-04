# bash completion for genomics
_genomics_completion()
{
    local cur prev words cword
    if declare -F _init_completion >/dev/null 2>&1; then
        _init_completion -n : || return
    else
        COMPREPLY=()
        cur="${COMP_WORDS[COMP_CWORD]}"
        prev="${COMP_WORDS[COMP_CWORD-1]}"
        words=("${COMP_WORDS[@]}")
        cword=$COMP_CWORD
    fi

    local commands="audit-configs audit-data convert snp-ancestry genomes-analyzer dataset-builders alphagenome genotype variant neural completion"
    local genotype="prepare-cache train evaluate pca-variance workbench sync-bcftools-artifacts single-gene-screen"
    local variant="materialize train evaluate analyze-counts"
    local convert="vcf-to-23andme"
    local snp="run"
    local genomes_analyzer="run"
    local dataset_builders="non-longevous"
    local non_longevous="build build-window visualize"
    local alphagenome="analyze integrate tracks"
    local neural="train test summarize pca-cache"
    local completion="bash"

    case "${words[1]}" in
        genotype) COMPREPLY=( $(compgen -W "$genotype" -- "$cur") ); return ;;
        variant) COMPREPLY=( $(compgen -W "$variant" -- "$cur") ); return ;;
        convert) COMPREPLY=( $(compgen -W "$convert" -- "$cur") ); return ;;
        snp-ancestry) COMPREPLY=( $(compgen -W "$snp" -- "$cur") ); return ;;
        genomes-analyzer) COMPREPLY=( $(compgen -W "$genomes_analyzer" -- "$cur") ); return ;;
        dataset-builders)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$dataset_builders" -- "$cur") ); return; fi
            if [[ ${words[2]} == non-longevous ]]; then COMPREPLY=( $(compgen -W "$non_longevous" -- "$cur") ); return; fi
            ;;
        alphagenome) COMPREPLY=( $(compgen -W "$alphagenome" -- "$cur") ); return ;;
        neural) COMPREPLY=( $(compgen -W "$neural" -- "$cur") ); return ;;
        completion) COMPREPLY=( $(compgen -W "$completion" -- "$cur") ); return ;;
    esac

    if [[ ${cword} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "$commands" -- "$cur") )
        return
    fi

    case "$prev" in
        --config|-c|config|--output|--json-output|--dataset-dir|--processed-dir|--results-dir|--checkpoint)
            if declare -F _filedir >/dev/null 2>&1; then
                _filedir
            else
                COMPREPLY=( $(compgen -f -- "$cur") )
            fi
            return
            ;;
    esac
}
complete -F _genomics_completion genomics
