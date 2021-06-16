#!/bin/bash

julia --project=. --startup-file=no --trace-compile=./../app_precompile.jl -e 'using Pkg; Pkg.test();'
julia --project=. -e 'using Pkg; Pkg.add("PackageCompiler"); using PackageCompiler; create_app(".", "./../CRISPRofftargetHunter_build"; precompile_statements_file = "./../app_precompile.jl", force = true);'

tar -czvf CRISPRofftargetHunter.tar.gz ./../CRISPRofftargetHunter_build