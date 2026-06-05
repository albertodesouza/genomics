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

    local commands="audit-configs audit-data config convert snp-ancestry genomes-analyzer dataset-builders alphagenome genotype variant neural completion"
    local config="describe schema validate"
    local genotype="prepare-cache split train test evaluate pca-variance workbench sync-bcftools-artifacts single-gene-screen"
    local variant="materialize train evaluate analyze-counts"
    local convert="vcf-to-23andme"
    local snp="run"
    local genomes_analyzer="run"
    local dataset_builders="non-longevous"
    local non_longevous="build build-window visualize"
    local alphagenome="analyze integrate tracks"
    local neural="train test summarize pca-cache"
    local completion="bash"

    _genomics_filedir()
    {
        if declare -F _filedir >/dev/null 2>&1; then
            _filedir "$@"
        else
            COMPREPLY=( $(compgen -f -- "$cur") )
        fi
    }

    _genomics_yaml_configs()
    {
        local item matches=()
        while IFS= read -r item; do
            matches+=("$item/")
        done < <(compgen -d -- "$cur")
        while IFS= read -r item; do
            case "$item" in
                *.yaml|*.yml) matches+=("$item") ;;
            esac
        done < <(compgen -f -- "$cur")
        COMPREPLY=("${matches[@]}")
    }

    if [[ ${cword} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "$commands" -- "$cur") )
        return
    fi

    case "${words[1]}" in
        genotype)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$genotype" -- "$cur") ); return; fi
            case "${words[2]}" in
                prepare-cache|split|train|test|evaluate|pca-variance|single-gene-screen) _genomics_yaml_configs; return ;;
                workbench|sync-bcftools-artifacts) _genomics_filedir; return ;;
            esac
            ;;
        variant)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$variant" -- "$cur") ); return; fi
            case "${words[2]}" in
                train|evaluate) _genomics_yaml_configs; return ;;
                materialize|analyze-counts) _genomics_filedir; return ;;
            esac
            ;;
        config)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$config" -- "$cur") ); return; fi
            case "${words[2]}" in
                describe|schema) COMPREPLY=( $(compgen -W "genotype variant" -- "$cur") ); return ;;
                validate) _genomics_yaml_configs; return ;;
            esac
            ;;
        convert)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$convert" -- "$cur") ); return; fi
            _genomics_filedir; return ;;
        snp-ancestry)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$snp" -- "$cur") ); return; fi
            _genomics_yaml_configs; return ;;
        genomes-analyzer)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$genomes_analyzer" -- "$cur") ); return; fi
            _genomics_yaml_configs; return ;;
        dataset-builders)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$dataset_builders" -- "$cur") ); return; fi
            if [[ ${words[2]} == non-longevous && ${cword} -eq 3 ]]; then COMPREPLY=( $(compgen -W "$non_longevous" -- "$cur") ); return; fi
            if [[ ${words[2]} == non-longevous ]]; then _genomics_yaml_configs; return; fi
            ;;
        alphagenome)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$alphagenome" -- "$cur") ); return; fi
            _genomics_filedir; return ;;
        neural)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$neural" -- "$cur") ); return; fi
            _genomics_yaml_configs; return ;;
        completion)
            if [[ ${cword} -eq 2 ]]; then COMPREPLY=( $(compgen -W "$completion" -- "$cur") ); return; fi
            ;;
    esac

    case "$prev" in
        --config|-c|config|validate|--output|--json-output|--dataset-dir|--processed-dir|--results-dir|--checkpoint)
            _genomics_filedir
            return
            ;;
    esac
}
complete -F _genomics_completion genomics
