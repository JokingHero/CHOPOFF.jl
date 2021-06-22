# CRISPRofftargetHunter

Software designed for efficient, precise and fast identification of all off-targets for given CRISPR guideRNA's. It leverages couple of alghoritms designed specifically for this purpose.

CRISPRofftargetHunter is **extensively tested**:

* unit tests for each function
* friction tests where three different implementations of the same functionality must report the same results
* end-to-end tests where we run whole pipeline on specially designed sample genome and compare results with CRISPRitz software

CRISPRofftargetHunter is designed specifically for CRISPR alignments of guideRNA's, **allowing for deletions, insretions and mismatches**. Implemented alghoritms allow you to find off-targets within **distance as large as you want**!

CRISPRofftargetHunter has an alghoritm (see: `build_sketchDB` and `search_sketchDB`) designed for **super fast estimation of number of off-targets** in the genome. These estimations can never report counts less than reality, but can only over-estimate! This can be used to quickly design libraries of guideRNA's for the entire genomes, as promising guides can be quickly sorted using CRISPRofftargetHunter.

CRISPRofftargetHunter has **support for multiple-cores**, we use standard julia configuration for that.

## Citation

TODO

## LICENSE

Copyright (C) 2021  Kornel Labun

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

## Self-contained build

It is possible to build CRISPRofftargetHunter into standalone application.
This can be achieved by running `./build_standalone.sh` script from the main directory. Script will
produce binary in a new folder outside the main directory. Then you can run from inside that folder `./bin/CRISPRofftargetHunter --help` or `./bin/CRISPRofftargetHunter build --help`.

## Main API

Run CRISPRofftargetHunter as an application. From the directory of the package run:

```bash
julia --threads 4 --project="." CRISPRofftargetHunter.jl --help
```

## Example

TODO

## Public Interface

You can also use CRISPRofftargetHunter as a normal julia package with the exported functions.

```@docs
Motif
build_linearDB
search_linearDB
build_sketchDB
search_sketchDB
build_treeDB
search_treeDB
inspect_treeDB
```
