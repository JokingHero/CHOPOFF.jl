[![Tests](https://github.com/JokingHero/CRISPRofftargetHunter.jl/actions/workflows/test.yml/badge.svg)](https://github.com/JokingHero/CRISPRofftargetHunter.jl/actions/workflows/test.yml)

# CRISPRofftargetHunter  


## About

Julia package that allows to search for guideRNA off-targets within genome
of interest while allowing for arbitrary distance
(allows for insertions, deletions and mismatches).

The goal is to make it fast and reliable.

## TODO

# Must do

* prepare pipeline for benchmark of speed against hg38v34!!! (~1 week)
* fix find_all function to not skip ending NNN (allow them) (1 d)
* code review and cleanup (3 d)
* VCF support (5 d)
* test for EVERY function (2 d)
* add end-to-end tests for bindDB and dictDB, friction tests too (1 d)
* more guides for tests - currently we have 20 - we want large range of guides with different properties (1 d)
* add tests for Cpf1 style PAMs (1 d)
* add optimizations: @inbounds, @inline etc. (2 d)
* add EARLY stopping!!! 
* more testing against CRISPRitz (2 d - 2 weeks)
* add interface for our standalone and main for binDB, dictDB etc. (1 d)
* publish code and package!!! - make sure it works
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
