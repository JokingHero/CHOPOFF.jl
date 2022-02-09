[![Tests & Release](https://github.com/JokingHero/CRISPRofftargetHunter.jl/actions/workflows/build_standalone.yml/badge.svg?branch=master)](https://github.com/JokingHero/CRISPRofftargetHunter.jl/actions/workflows/build_standalone.yml)

# CRISPRofftargetHunter  


## About

Julia package that allows to search for guideRNA off-targets within genome
of interest while allowing for arbitrary distance
(allows for insertions, deletions and mismatches).

The goal is to make it fast and reliable.

## TODO

# Must do

* fmi index searches
    * lossless seed  + indexed pam loci -> search_pamDB
    * partially seed on the prefix and then + lossless seed?/alignment -> search_fmiDB_raw
* improve motifDB with chunked lossless seed?
* what is the probability that guide and off-target will be within certain distance, based only on some sort of transform of their letters - or rather figure out 0 probability in a fast way!


* test small motifDB with UInt128, test cashe on hgv38!, test search_pamDB

* add EARLY stopping!!! 
* VCF support (5 d) 
* Loci Range alternative where we have default value of just one UInt32, else we have UInt8/UInt16 as width! - will decrease the database sizes considerably because most of the guides are unique!

* bloom filter paralelization?!
* implement xor or ribbon filters with upgrade that allows 
* prepare pipeline for benchmark of speed against hg38v34!!! (~1 week)
* code review and cleanup (3 d)
* test for EVERY function (2 d)
* more guides for tests - currently we have 20 - we want large range of guides with different properties (1 d)
* add tests for Cpf1 style PAMs (1 d)
* add optimizations: @inbounds, @inline etc. (2 d)
* more testing against CRISPRitz (2 d - 2 weeks)
* publish code and package!
* finish writing paper!!!

# Maybe

* paralelize binDB! sketchDB (aromic arrays of some sort)?
* implement faster, non-recurent version of comb_of_d


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
