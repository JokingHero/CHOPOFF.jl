#!/usr/bin/env bash
# Source this file to configure the copied development installation.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Source this script: source scripts/dev_env.sh" >&2
    exit 1
fi

_chopoff_setup_env() {
    local repo_root julia_bin julia_dir
    repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || return 1
    julia_bin="${JULIA_BIN:-$repo_root/../Soft/bin/julia}"
    if [[ ! -x "$julia_bin" ]]; then
        echo "Julia executable not found: $julia_bin. Set JULIA_BIN to a full Julia installation." >&2
        return 1
    fi
    julia_dir="$(cd "$(dirname "$julia_bin")" && pwd)" || return 1
    if [[ -z "${JULIA_DEPOT_PATH:-}" && ! -d "$repo_root/../Soft/julia_depot" ]]; then
        echo "Copied depot not found. Set JULIA_DEPOT_PATH to your package depot." >&2
        return 1
    fi
    if [[ ! -f "${CHOPOFF_PROFILE_ENV:-$repo_root/../profiletools}/Project.toml" ]]; then
        echo "Profiling environment not found. Set CHOPOFF_PROFILE_ENV to its directory." >&2
        return 1
    fi
    export JULIA_BIN="$julia_dir/$(basename "$julia_bin")"
    case ":$PATH:" in
        *":$julia_dir:"*) ;;
        *) export PATH="$julia_dir:$PATH" ;;
    esac
    export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$repo_root/../Soft/julia_depot:}"
    export CHOPOFF_PROFILE_ENV="${CHOPOFF_PROFILE_ENV:-$repo_root/../profiletools}"
}

if _chopoff_setup_env; then
    unset -f _chopoff_setup_env
else
    unset -f _chopoff_setup_env
    return 1
fi
