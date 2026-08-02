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

Building of `prefixHashDB` database for standard Cas9 `--motif` with support for up to levenshtein distance 3 `--distance` for an example genome using 10 threads.

```bash
export JULIA_NUM_THREADS=10  
EXAMPLE_GENOME="./test/sample_data/genome/semirandom.fa"
CHOPOFF build --name Cas9_hg38 --genome "$EXAMPLE_GENOME" -o out_dir/phDB_16_3/ --distance 3 --motif Cas9 prefixHashDB
```

Searching of above database for all off-targets for guides listed in `--guides` up to the 2 levenshtein distance `--distance` using 15 threads, writing the results into `--output` file. Because `--early_stopping` argument is not supplied below, by default
`prefixHashDB` will search for up to 1e6 off-targets per guide per distance. Pay attention that default guides for the Cas9, are 20bp long, as can be inspected in the example file.

```bash
export JULIA_NUM_THREADS=15  
EXAMPLE_GUIDES="./test/sample_data/guides.txt"
CHOPOFF search --database phDB_16_3/ --guides "$EXAMPLE_GUIDES" --output out_dir/phDB_16_2.csv --distance 2 prefixHashDB
```

## No-build application

Run CHOPOFF package as an **application**, without building first. From the directory of the package run:

```bash
julia --threads 4 --project="." ./src/CHOPOFF.jl --help
```

## Quick Use

For search of off-targets you have a couple of options:
- **[prefixHashDB](@ref "search_prefixHashDB")** - the fastest, we apply hashes to symbolic alignments for fast filtering of OTs  
- [prefixHashScan](prefix_hash_scan.md) - direct indexed-FASTA or 2bit search for registered or custom motifs, including PAMless searches
- [linearDB](@ref "search_linearDB") - most rigorously tested
- [motifDB](@ref "search_motifDB") - on top of linearDB we apply pigeonhole principle like filter which you can adjust
- [treeDB](@ref "search_treeDB") - will work best for longer gRNAs, uses vantage point for filteirng
- [Sassy search](sassy_search.md) - direct genome search with SIMD minima and DP traceback (no prebuilt CHOPOFF DB)

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

Copyright © 2025 Kornel Labun

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
