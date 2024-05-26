[![Tests & Release](https://github.com/JokingHero/CHOPOFF.jl/actions/workflows/build_standalone.yml/badge.svg?branch=master)](https://github.com/JokingHero/CHOPOFF.jl/releases/tag/latest) 
![Coverage is High](./coverage/coverage_fraction.svg) 
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jokinghero.github.io/CHOPOFF.jl/) 

# 

![Logo](./docs/src/assets/logo-dark.png#gh-dark-mode-only)
![Logo](./docs/src/assets/logo.png#gh-light-mode-only)
## <p align="center">Sensitive and fast CRISPR off-target detection</p>

## About CHOPOFF

* many fast alignment algorithms optimized specifically for CRISPR
* search for larger distances allowing for mismatches and bulges
* support for ambiguous bases
* arbitrarily large genomes
* VCF support - with multiple overlapping SNPs
* off-target filtering that is near-instant and alignment-free 
* pruning of off-targets by their location (remove overlapping, competing off-targets)
* extensively tested
* full framework that can be extended for your own algorithms with ease


## Requirements

* Some algorithms generate as many files as there are prefixes (e.g. for prefix 7 - this will make 4^7 - 16384 files). This strategy allows us to operate the searches independently on multiple cores and not get throttled when querying large number of the guides. However, some systems have artificial limits on the number of open files, for example in Ubuntu 'ulimit -n' will show the limit. Increase the limits, if it creates problems for you.

* When using many cores for building the indexes - you have to have around ~1 GB of RAM per thread.


## Install and use as a standalone application

It is possible to build CHOPOFF into standalone application - which includes all dependencies and Julia into one compiled software. This is **recommended** method for using of CHOPOFF when you are not a developer. If you know how to code in Julia, you might make use of the whole framework using CHOPOFF as a package (see docs [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jokinghero.github.io/CHOPOFF.jl/)).

To build a standalone application run `./build_standalone.sh` script from the main directory. Script will run all tests and also 
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

## R integration with crisprVerse

Visit [crisprCHOPOFF](https://github.com/JokingHero/crisprCHOPOFF) - R package that allows you to use CHOPOFF with crisprVerse.  


## Support

You can buy me a [coffee](https://www.buymeacoffee.com/kornellabun) to show some love and appreciation!

<img src="./docs/src/assets/bmc_qr.png" width="25%"/>


## License  

Copyright © 2024 Kornel Labun


The software is provided under the following terms:

1. Academic and Non-Commercial Use License: The software is licensed under the aGPL-3.0 License, attached here inside Non_Commercial_LICENSE file. If not, see <https://www.gnu.org/licenses/>. This license applies only to academic institutions and for non-commercial use. Users wishing to use the software for commercial purposes must obtain a separate license.

2. Commercial Use License: Parties interested in using the software for commercial purposes are invited to contact Vestlandets Innovasjonsselskap AS (hei@visinnovasjon.no) for a commercial use license.