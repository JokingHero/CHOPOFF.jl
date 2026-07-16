```@meta
CollapsedDocStrings = true
```

# prefixHashScan search

`search_prefixHashScan` searches an indexed FASTA reference directly without
building a CHOPOFF genome database. Separate Cas9 and Cas12a distance-3 kernels
amortize one genome scan across all supplied guides.

## Supported configuration

- 1-64 unambiguous guides;
- Cas9: 20-base guide with an `NGG` PAM;
- Cas12a: 21-base guide with a `TTTV` PAM;
- edit distance 3;
- 16-base symbolic prefix filter;
- FASTA reference with a standard `.fai` index;
- strict ACGT candidate windows.

Any candidate whose complete 23-base Cas9 or 25-base Cas12a guide/PAM window
contains a non-ACGT base is skipped. Ambiguous query guides are rejected.

## API

```@docs
search_prefixHashScan
```

```julia
using CHOPOFF, BioSequences

guides = LongDNA{4}.([
    "ACGTACGTACGTACGTACGT",
    "TGCATGCATGCATGCATGCA",
])

search_prefixHashScan(
    guides,
    "genome.fa",
    "offtargets.csv";
    scan_threads = Threads.nthreads(),
    verbose = true,
)

# Cas12a
search_prefixHashScan(
    LongDNA{4}.(["TGCATGCATGCATGCATGCAT"]),
    "genome.fa",
    "cas12a_offtargets.csv";
    motif = "Cas12a",
)
```

The output uses the normal CHOPOFF detail columns:
`guide`, `alignment_guide`, `alignment_reference`, `distance`, `chromosome`,
`start`, and `strand`.

With `verbose=true`, the search reports the resolved scan backend, lookup mode,
query-build mode, scheduler, thread count, and chunk size. Reporting does not
enable hot-loop statistics.

## Command-line usage

```bash
CHOPOFF search --distance 3 \
  --guides guides.txt \
  --output offtargets.csv \
  prefixHashScan --genome genome.fa --motif Cas12a
```

The CLI always reports the resolved execution mode. It uses Julia's configured
thread count, for example `JULIA_NUM_THREADS=12`.

## Implementation boundary

The three-argument method defaults to Cas9 and accepts `motif="Cas12a"`. The
four-argument `Motif` method supports both standard geometries; its tuning
keywords and generic legacy fallback remain experimental. Geometry dispatch
happens before separate Cas9 and Cas12a hot loops.
