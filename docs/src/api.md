# Public Interface

When using CRISPRofftargetHunter as a normal julia package with the exported functions we provide many methods for working with CRISPR gRNA and off-targets for your convenience.

## Abstract gRNA

To encapsulate motifs that we search for on the genome, we can use abstract definition of gRNA/off-target, which is defined in `Motif` data structure.

```@docs
Motif
length_noPAM
setdist
setambig
```

## Find off-targets

```@docs
DBInfo
gatherofftargets!
```

## Align gRNA and off-target

```@docs
isinclusive
hamming
levenshtein
Aln
align
```

## Utils

### gRNAs and kmers
```@docs
getseq
as_kmers
as_skipkmers
all_kmers
minkmersize
``` 

### Persistence
```@docs
save
load
```

## Alignment-free filters for gRNAs

```@docs
build_hashDB
search_hashDB
```

## Find all off-targets

### Prefix-Suffix partial alignment

```@docs
build_linearDB
search_linearDB

build_motifDB
search_motifDB
```

### Vantage-Point tree

```@docs
build_treeDB
inspect_treeDB
search_treeDB
```