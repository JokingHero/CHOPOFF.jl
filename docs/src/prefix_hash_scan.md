```@meta
CollapsedDocStrings = true
```

# prefixHashScan search

`search_prefixHashScan` searches an indexed FASTA reference directly without
building a CHOPOFF genome database. The supported entrypoint is specialized for
Cas9 at edit distance 3 and amortizes one genome scan across all supplied
guides.

## Supported configuration

- 1-64 unambiguous 20-base Cas9 guides;
- Cas9 `NGG` PAM and edit distance 3;
- 16-base symbolic prefix filter;
- FASTA reference with a standard `.fai` index;
- strict ACGT candidate windows.

Any candidate whose complete 23-base guide/PAM window contains a non-ACGT base
is skipped. Ambiguous query guides are rejected.

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
  prefixHashScan --genome genome.fa
```

The CLI always reports the resolved execution mode. It uses Julia's configured
thread count, for example `JULIA_NUM_THREADS=12`.

## Implementation boundary

The documented three-argument method above is the supported interface. The
four-argument method accepting a `Motif` and tuning keywords remains an
experimental benchmark/parity engine. Cas12a and other geometries should use
separate specialized kernels rather than adding motif branches to the Cas9 hot
loop.
