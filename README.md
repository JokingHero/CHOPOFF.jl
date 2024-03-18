[![Tests & Release](https://github.com/JokingHero/CHOPOFF.jl/actions/workflows/build_standalone.yml/badge.svg?branch=master)](https://github.com/JokingHero/CHOPOFF.jl/releases/tag/latest) 
![Coverage is High](./coverage/coverage_fraction.svg) 
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jokinghero.github.io/CHOPOFF.jl/) 

# CHOPOFF

![Logo](./docs/src/assets/logo-dark.png#gh-dark-mode-only)
![Logo](./docs/src/assets/logo.png#gh-light-mode-only)

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


## Support

You can buy me a [coffee](https://www.buymeacoffee.com/kornellabun) to show some love and appreciation!

<img src="./docs/src/assets/bmc_qr.png" width="25%"/>


## License  

Copyright © 2022 Kornel Labun

License for non-commercial applications is aGPL-3.0. 
For commercial applications you should acquire permission or licensing contract.

https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0)

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
