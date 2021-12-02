#!/bin/bash

julia --project=. --startup-file=no --trace-compile=./precompile/app_precompile.jl ./precompile/precompile.jl
julia --project=. -e 'using Pkg; Pkg.add("PackageCompiler"); using PackageCompiler; create_app(".", "./build/"; precompile_statements_file = "./precompile/app_precompile.jl", force = true);'
