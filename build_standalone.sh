#!/bin/bash

if [ -z "$1" ]
  then
    M_T=`uname -m`
    S_T=`uname -s`
    C_V=`julia --project=. -e 'using PkgVersion; print(string(@PkgVersion.Version));'`
    J_V=`julia --version`
    J_V=($J_V)
    N="CRISPRofftargetHunter_v${C_V}_julia_v${J_V[-1]}_${S_T}_${M_T}"
  else
    N=${1}
fi

julia --project=. --startup-file=no --trace-compile=./../app_precompile.jl -e 'using Pkg; Pkg.test();'
julia --project=. -e 'using Pkg; Pkg.add("PackageCompiler"); using PackageCompiler; create_app(".", "./../'${N}'"; precompile_statements_file = "./../app_precompile.jl", force = true);'

tar -czvf "./../${N}.tar.gz" "./../${N}"
