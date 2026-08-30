```@meta
CollapsedDocStrings = true
```

# Find all off-targets

## Fastest method with symbolic alignments and prefix hashes

```@docs
build_prefixHashDB
search_prefixHashDB
```

## Prefix-Suffix partial alignment

```@docs
build_linearDB
search_linearDB

build_motifDB
search_motifDB
```

## Vantage-Point tree

```@docs
build_treeDB
inspect_treeDB
search_treeDB
```

## Path templates

```@docs
build_PathTemplates
```

## prefixHashScan direct genome search

prefixHashScan is a direct search mode (no prebuilt CHOPOFF DB needed) that
reuses prefixHashDB's symbolic prefix paths and scans an indexed FASTA or 2bit
reference with SIMD kernels. It handles registered motifs, custom motifs and
PAMless searches. See [prefixHashScan search](prefix_hash_scan.md) for the
supported configurations, CLI usage and tuning options.

```@docs
search_prefixHashScan
```

## Sassy direct genome search

Sassy is a direct search mode (no prebuilt CHOPOFF DB needed) that combines SIMD
candidate generation with DP traceback for full alignments. See
[Sassy search](sassy_search.md) for algorithm and usage details.

## VCF

```@docs
build_vcfDB
search_vcfDB
```
