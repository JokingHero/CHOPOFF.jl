```@meta
CollapsedDocStrings = true
```

# prefixHashScan search

`search_prefixHashScan` searches an indexed FASTA or 2bit reference directly
without building a CHOPOFF genome database. Hand-written Cas9 and Cas12a
kernels remain the canonical fast paths; a motif-specialized generic kernel
handles other PAM sequences, PAM positions, strand selections, and motifs
without a PAM.

## Supported configuration

- any number of unambiguous guides, processed in batches of at most 64;
- any registered motif name or custom `Motif` object;
- one contiguous PAM block at any position, or no PAM;
- edit distance 0, 1, 2, 3, or 4;
- 0 through 3 ambiguous IUPAC reference bases per guide/PAM window;
- a symbolic prefix of up to 16 bases;
- FASTA reference with a standard `.fai` index, or a `.2bit` reference without
  a sidecar index;

The typed generic SIMD engine is selected when the prefix is 16 bases, the
guide is 16 through 64 bases, the complete motif spans at most 65 bases, and
the guide retains at least 16 bases at the requested distance
(`guide length - distance >= 16`). Other valid motifs use the exact legacy
backend. Guide lists may be larger than 64 because the public method batches
them automatically.

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
    "genome.2bit",
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

# No PAM: an all-X PAM description produces an empty PAM range.
pamless = Motif(
    "pamless", repeat("N", 20), repeat("X", 20),
    true, true, 3, true, 0)
search_prefixHashScan(
    guides,
    "genome.fa",
    "pamless_offtargets.csv";
    motif = pamless,
    distance = 3,
)

# Variable-length guide: 25 guide bases followed by an NNT PAM.
guide25_nnt = Motif(
    "25N_NNT",
    repeat("N", 25) * "NNT",
    repeat("X", 25) * "NNT",
    true, true, 4, true, 0,
)
search_prefixHashScan(
    LongDNA{4}.(["ACGTACGTACGTACGTACGTACGTA"]),
    "genome.fa",
    "guide25_nnt_offtargets.csv";
    motif = guide25_nnt,
    distance = 3,
)
```

The 25-base example resolves to the typed generic fast engine, not `:legacy`.
Only 20- and 21-base guides currently have shipped symbolic path assets. Other
eligible lengths generate paths once before guide batching and then use the
same optimized scanner.

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

FASTA preserves all listed IUPAC symbols. The 2bit format preserves ambiguity
only as `N` blocks, so other symbols converted to 2bit are handled as `N`.
Soft-mask blocks do not suppress candidates; prefixHashScan searches them as
their uppercase bases.

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
  prefixHashScan --genome genome.2bit --motif Cas12a --ambig_max 3
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

## Generic engine overhead

A local single-thread microbenchmark forced the identical Cas9/distance-3 motif
through the hand-written and typed generic scanners over 32 MB of deterministic
random A/C/G/T reference. Eleven alternating timed runs gave:

| Scanner | Median | Throughput |
|---|---:|---:|
| Hand-written Cas9 | 66.7 ms | 480 MB/s |
| Forced typed generic | 78.6 ms | 407 MB/s |

Candidate counts were identical. The generic scanner took 17.8% longer and
retained about 85% of specialized Cas9 throughput. This isolates motif scanning,
prefix packing, and failed hash lookup; end-to-end overhead varies with PAM
frequency, query hits, verification, I/O, and output volume.

## Implementation boundary

The three-argument method defaults to Cas9/distance 3 and accepts any registered
motif name or custom `Motif` plus `distance=0:4`. Geometry dispatch happens
before the canonical or typed generic hot loop. Symbolic path assets are reused
by guide length, distance, and prefix length; scanner specialization remains
independent of path generation.
