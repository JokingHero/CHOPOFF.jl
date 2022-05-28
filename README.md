[![Tests & Release](https://github.com/JokingHero/CRISPRofftargetHunter.jl/actions/workflows/build_standalone.yml/badge.svg?branch=master)](https://github.com/JokingHero/CRISPRofftargetHunter.jl/releases/tag/latest) 
![Coverage is High](./coverage/coverage_fraction.svg)

# CRISPRofftargetHunter  

## About

Julia framework related to CRISPR off-targets:
* fast alignment alghoritms optimized for CRISPR
* find quickly all off-targets on the arbitrairly large genomes, with arbitrary distances, and ambiogous bases support
* VCF support - with multiple overlapping SNPs
* near-instant alignment-free off-target filtering
* pruning of off-targets by their location (remove overlapping, competing off-targets)
* extensively tested

We have chosen Julia because it is possible to quickly prototype new
alghoritms and not compromise the speed, but also, Julia's syntax is python like and it is easy to jump right into coding.

Documentation is available here.


## Requirements

* Some of the alghoritms generate as many files as there are prefixes (e.g. for prefix 7 - this will make 4^7 - 16384) when using '--detail' option. This strategy allows us to operate the searches independently on multiple cores and not get throtled when querying large number of the guides. However, some systems have artificial limits on the number of open files, for example in ubuntu 'ulimit -n' will show the limit. Increase the limits, if it creates problems for you.

* When using many cores for building the indexes - you have to have around ~1GB of RAM per thread.


## License  

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
