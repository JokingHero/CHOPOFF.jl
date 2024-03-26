# CHOPOFF.jl

## About

Uncompromising finding of CRISPR off-targets:
* many fast alignment algorithms optimized specifically for CRISPR
* search for larger distances allowing for mismatches and bulges
* support for ambiguous bases
* arbitrarily large genomes
* VCF support - with multiple overlapping SNPs
* near-instant alignment-free off-target filtering
* pruning of off-targets by their location (remove overlapping, competing off-targets)
* extensively tested
* full framework that can be extended for your own algorithms with ease

## Requirements

* Some algorithms generate as many files as there are prefixes (e.g. for prefix 7 - this will make 4^7 - 16384 files). This strategy allows us to operate the searches independently on multiple cores and not get throttled when querying large number of the guides. However, some systems have artificial limits on the number of open files, for example in Ubuntu 'ulimit -n' will show the limit. Increase the limits, if it creates problems for you.

* When using many cores for building the indexes - you have to have around ~1 GB of RAM per thread.

## Standalone application

It is possible to build CHOPOFF into standalone application - which includes all dependencies and Julia into one compiled software. This is **recommended** method for using of CHOPOFF when you are not a developer. If you know how to code in Julia, you might make use of the whole framework using CHOPOFF as a package.

To build a standalone application run `./build_standalone.sh` script from the main directory. Script will
produce binary in a "build" folder. Then you can run from inside that folder `./bin/CHOPOFF --help`. To learn about building a database run `./bin/CHOPOFF build --help` and to use existing database check out `./bin/CHOPOFF search --help`. It is possible to skip testing + precompile step to speed up the build process with `./build_standalone.sh --noprecompile`.

You can alternatively download the latest release from the releases' page on the GitHub.

When using application as self-contained compiled software, you can control number of cores by setting `JULIA_NUM_THREADS` environment variable.

**Example commands for using standalone**

Building of `prefixHashDB` database for standard Cas9 `--motif` with support for up to levenshtein distance 4 `--distance` for an example genome using 10 threads.

```bash
export JULIA_NUM_THREADS=10  
CHOPOFF build --name Cas9_hg38 --genome hg38.fa -o out_dir/phDB_16_4/ --distance 4 --motif Cas9 prefixHashDB --hash_length 16
```

Searching of above database for all off-targets for guides listed in `--guides` up to the 2 levenshtein distance `--distance` using 15 threads, writing the results into `--output` file. Because `--early_stopping` argument is not supplied below, by default
`prefixHashDB` will search for up to 1e6 off-targets per guide per distance.

```bash
export JULIA_NUM_THREADS=15  
CHOPOFF search --database phDB_16_4/ --guides guides_wo_PAM.txt --output out_dir/phDB_16_2.csv --distance 2 prefixHashDB
```

## No-build application

Run CHOPOFF package as an **application**, without building first. From the directory of the package run:

```bash
julia --threads 4 --project="." ./src/CHOPOFF.jl --help
```

## Quick Use

For search of off-targets you have a couple of options:
- **[prefixHashDB](@ref "search_prefixHashDB")** - the fastest, we apply hashes to symbolic alignments for fast filtering of OTs  
- [linearDB](@ref "search_linearDB") - most rigorously tested
- [motifDB](@ref "search_motifDB") - on top of linearDB we apply pigeonhole principle like filter which you can adjust
- [treeDB](@ref "search_treeDB") - will work best for longer gRNAs, uses vantage point for filteirng
- [fmiDB](@ref "search_fmiDB") - the smallest file size, very fast, but only for distances ≤ 2
- [binaryFuseFilterDB](@ref "search_binaryFuseFilterDB") - uses hashing on top of FM-index

If you would like to **filter** or **rank** gRNAs to only those that are most likely off-target free you want to use [hashDB](@ref "search_hashDB") or
slower and much larger, but less probabilistic [dictDB](@ref "search_dictDB").

For VCF file support use [vcfDB](@ref "search_vcfDB"). 
For use of the framework as a Julia package consult the documentation - Public Interface section.


## Support

You can buy me a [coffee](https://www.buymeacoffee.com/kornellabun) to show some love and appreciation!

```@raw html
<img src="./assets/bmc_qr.png" width="25%"/>
```


## LICENSE

Copyright © 2022 Kornel Labun

License for non-commercial applications is aGPL-3.0. 
For commercial applications you should acquire permission or licensing contract.

<https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0)>

This program is free software for non-commercial applications: 
you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program inside LICENSE file. 
If not, see <https://www.gnu.org/licenses/>.
