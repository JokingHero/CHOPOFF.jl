# ARTEMIS

## About

Uncompromised finding of CRISPR off-targets:
* many fast alignment alghoritms optimized specifically for CRISPR
* search for larger distances allowing for mismatches and bulges
* support for ambiguous bases
* arbitrairly large genomes
* VCF support - with multiple overlapping SNPs
* near-instant alignment-free off-target filtering
* pruning of off-targets by their location (remove overlapping, competing off-targets)
* extensively tested
* full framework that can be extended for your own alghoritms with ease

## Requirements

* Some of the alghoritms generate as many files as there are prefixes (e.g. for prefix 7 - this will make 4^7 - 16384 files). This strategy allows us to operate the searches independently on multiple cores and not get throtled when querying large number of the guides. However, some systems have artificial limits on the number of open files, for example in ubuntu 'ulimit -n' will show the limit. Increase the limits, if it creates problems for you.

* When using many cores for building the indexes - you have to have around ~1GB of RAM per thread.

## Build application

It is possible to build ARTEMIS into standalone application - which includes all dependencies and Julia into one compiled software. This is **recommended** method for using of ARTEMIS when you are not a developer. If you know how to code in Julia, you might make use of the whole framework using ARTEMIS as a package.

To build a standalone application run `./build_standalone.sh` script from the main directory. Script will
produce binary in a "build" folder. Then you can run from inside that folder `./bin/ARTEMIS --help`. To learn about building a database run `./bin/ARTEMIS build --help` and to use existing database check out `./bin/ARTEMIS search --help`. It is possible to skip testing + precompile step to speed up the build process with `./build_standalone.sh --noprecompile`.

You can alternatively download latest release from the releases page on the github.

When using application as self-contained compiled software, you can control number of cores by setting `JULIA_NUM_THREADS` environment variable.

## No-build application

Run ARTEMIS package as an **application**,  without building first. From the directory of the package run:

```bash
julia --threads 4 --project="." ./src/ARTEMIS.jl --help
```

## Quick Use

For search of off-targets you have a couple of options:
- `linearDB` - most rigorously tested
- `motifDB` - the fastest
- `treeDB` - will work best for longer gRNAs
- `fmiDB` - smallest file size, very fast, but only for distances ≤ 2

If you would like to **filter** or **rank** gRNAs to only those that are most likely off-target free you want to use `hashDB` or
slower and much larger, but less probabilistic `dictDB`.

For VCF file support use `vcfDB`. 
For use of the framework as a Julia package consult the documentation - Public Interface section.

## LICENSE

Copyright (C) 2022  Kornel Labun

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
