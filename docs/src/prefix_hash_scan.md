```@meta
CollapsedDocStrings = true
```

# prefixHashScan search

`search_prefixHashScan` searches an indexed FASTA reference directly without
building a CHOPOFF genome database. Separate Cas9 and Cas12a kernels serve edit
distances 0 through 4 and amortize one genome scan across all supplied guides.

## Supported configuration

- 1-64 unambiguous guides;
- Cas9: 20-base guide with an `NGG` PAM;
- Cas12a: 21-base guide with a `TTTV` PAM;
- edit distance 0, 1, 2, 3, or 4;
- 0 through 3 ambiguous IUPAC reference bases per guide/PAM window;
- 16-base symbolic prefix filter;
- FASTA reference with a standard `.fai` index;
- any number of guides, processed in batches of at most 64.

Ambiguous query guides are rejected. `ambig_max=0` remains the default and
skips candidate windows containing ambiguous reference bases.

Distance 4 is the supported upper-bound and benchmarking mode, not a tuned
fast path. Its p16 query is substantially larger than d0-d3; a 61-guide Cas9
query on the development host took 10.8 seconds to build, occupied 1.09 GB, and
reached 4.51 GB peak process RSS.

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
    distance = 2,
    scan_threads = Threads.nthreads(),
    verbose = true,
)

# Cas12a
search_prefixHashScan(
    LongDNA{4}.(["TGCATGCATGCATGCATGCAT"]),
    "genome.fa",
    "cas12a_offtargets.csv";
    motif = "Cas12a",
    distance = 1,
)

# Allow up to three ambiguous reference bases in each guide/PAM window.
ambiguous_motif = Motif("Cas9"; distance = 3, ambig_max = 3)
search_prefixHashScan(
    guides,
    "genome.fa",
    "ambiguous_offtargets.csv";
    motif = ambiguous_motif,
    distance = 3,
)
```

The output uses the normal CHOPOFF detail columns:
`guide`, `alignment_guide`, `alignment_reference`, `distance`, `chromosome`,
`start`, and `strand`.

With `verbose=true`, the search reports the resolved scan backend, lookup mode,
query-build mode, scheduler, thread count, and chunk size. Reporting does not
enable hot-loop statistics.

## Reference ambiguity semantics

`ambig_max` limits ambiguous reference positions across the complete candidate
window: 23 bases for Cas9 or 25 bases for Cas12a. Supported IUPAC symbols are
`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, and `N`, including lowercase
input. Each ambiguous position counts once, regardless of how many nucleotides
its code represents.

- A window with more than `ambig_max` ambiguous positions is skipped.
- PAM ambiguity must be compatible with the motif. For example, `ARG` is
  compatible with the Cas9 `NGG` PAM, while `AYG` is not.
- A compatible ambiguous guide base matches at zero edit cost. Incompatible
  possibilities are handled by the normal edit-distance calculation.
- The ambiguity allowance is separate from edit distance: distance 3 and
  `ambig_max=3` permit both limits in the same candidate.
- Unsupported reference symbols such as `X` are skipped.

For an ambiguous 16-base prefix, prefixHashScan enumerates only its compatible
two-bit hash variants, then merges their guide masks. The complete candidate is
materialized once and matching guide pairs are verified with IUPAC-aware Myers
distance and traceback. With at most three ambiguous positions, prefix
expansion is bounded by 64 variants for `N,N,N`; ambiguities outside the prefix
do not expand the hash lookup.

## Command-line usage

```bash
CHOPOFF search --distance 4 \
  --guides guides.txt \
  --output offtargets.csv \
  prefixHashScan --genome genome.fa --motif Cas12a --ambig_max 3
```

The CLI always reports the resolved execution mode. It uses Julia's configured
thread count, for example `JULIA_NUM_THREADS=12`.

## Cas9 ambiguity benchmark

The following development-host measurement used GRCh38, the 61 standard Cas9
guides in `test/local_human/data/guides_for_tests.txt`, distance 3, eight Julia
threads, two warmups, and 11 interleaved timed repetitions. Times are medians;
relative overhead is more portable than the absolute seconds.

| `ambig_max` | Prepared scan | Scan overhead | End to end | End-to-end overhead | Rows |
|------------:|--------------:|--------------:|-----------:|--------------------:|-----:|
| 0 | 1.922 s | control | 2.317 s | control | 25,826 |
| 1 | 2.259 s | +17.5% | 2.658 s | +14.7% | 25,827 |
| 2 | 2.271 s | +18.2% | 2.684 s | +15.9% | 25,827 |
| 3 | 2.257 s | +17.4% | 2.637 s | +13.8% | 25,827 |

All lower-ambiguity result multisets were contained in the next level. The one
additional row was already enabled at `ambig_max=1`; levels 2 and 3 added no
rows for this reference and guide set. Maximum end-to-end coefficient of
variation was 4.55%. Median allocation increased from 870.6 MB at level 0 to
878.0 MB at level 3, approximately 0.8%.

Reproduce or adapt the run with:

```bash
CHOPOFF_TUNING_STAGE=final \
CHOPOFF_TUNING_AMBIG_MAXES=0,1,2,3 \
CHOPOFF_TUNING_THREADS=8 \
CHOPOFF_TUNING_WARMUPS=2 \
CHOPOFF_TUNING_RUNS=11 \
CHOPOFF_TUNING_ALLOCATION_RUNS=3 \
CHOPOFF_TUNING_OUT=/tmp/prefix_hash_scan_ambiguity \
JULIA_NUM_THREADS=8 \
julia --project=. scripts/benchmark_prefix_hash_scan_tuning.jl
```

The runner writes raw timings, allocations, internal diagnostics, a summary
with control-relative ratios, final result files, and a multiset-inclusion
parity report.

## Implementation boundary

The three-argument method defaults to Cas9/distance 3 and accepts
`motif="Cas12a"` plus `distance=0:4`. The four-argument `Motif` method supports
both standard geometries; its tuning keywords and generic legacy fallback
remain experimental. Geometry dispatch happens before separate Cas9 and
Cas12a hot loops.
