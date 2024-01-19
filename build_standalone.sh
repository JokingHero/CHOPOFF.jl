#!/bin/bash

noprecompile=0
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -nop|--noprecompile) noprecompile=1 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ $noprecompile == 0 ]]
then
    julia --project=. --startup-file=no --trace-compile=./precompile/app_precompile.jl ./precompile/precompile.jl
fi
julia --project=. -e 'using Pkg; Pkg.add("PackageCompiler"); using PackageCompiler; create_app(".", "./build/"; precompile_statements_file = "./precompile/app_precompile.jl", force = true, include_lazy_artifacts = true, incremental = false, filter_stdlibs=true,);'
