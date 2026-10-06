# prefixHashScan

## Status and scope

`prefixHashScan` is an indexless CRISPR off-target search. Its public Julia API
accepts registered motifs and custom `Motif` objects at edit distances 0 through
4. Every eligible motif, Cas9 and Cas12a included, runs on one
motif-specialized generic kernel. Configurations outside that kernel's envelope
use the exact legacy engine. The search reuses the symbolic prefix paths from
`prefixHashDB`, builds a guide-specific query structure in memory, and scans the
reference genome directly.

"Indexless" means that no CHOPOFF genome database must be built. The FASTA
reader still requires the small, standard `.fai` random-access index. `.2bit`
references need no sidecar index.

The optimized envelope is:

- any motif with one contiguous PAM block at any position, or no PAM;
- 16 through 64 guide bases, complete motif span at most 65 bases;
- forward, reverse, or both strands;
- edit distance 0 through 4;
- a fixed 16-base prefix, with `guide length - distance >= 16`;
- 64 guides per query batch; larger lists are batched automatically;
- FASTA with a standard `.fai`, or a `.2bit` reference;
- reference windows with zero through three IUPAC-ambiguous positions
  (`motif.ambig_max`); query guides must be unambiguous;
- any CPU: x86 CPUs with AVX2 and BMI2 use `:avx2`, qualified CPUs use
  AVX-512F/BW, and every other CPU uses the `:portable` backend.

Production defaults:

- 2 MiB globally scheduled chunks;
- a 26-bit presence prefilter and 11-base directory bucket, fixed as
  `PREFIX_HASH_SCAN_PREFILTER_BITS` and `PREFIX_HASH_SCAN_BUCKET_BASES`;
- radix-ordered compact-directory lookup for prefilter survivors;
- per-guide hash construction on at most `scan_threads` tasks, followed by a
  deterministic serial heap merge;
- buffered scan, then verify, per chunk.

Two output modes exist. `detail` writes one aligned row per off-target.
`counts` writes `guide,D0,...,Dk,complete` per unique guide and skips traceback.
Guide lists larger than 64 run as sequential 64-guide batches. Detail mode
appends batches to a sibling staging file and atomically renames it after every
batch succeeds. Count mode merges batch matrices and writes one row per unique
guide after applying deterministic per-distance caps.

### Source layout

- `src/db_prefix_hash_scan.jl`: shared types, plan resolution, orchestration,
  legacy engine, and public API;
- `src/prefix_hash_scan/query.jl`: symbolic paths, hashes, directory, prefilter;
- `src/prefix_hash_scan/isa.jl`: CPU feature detection, backend resolution, and
  the backend-specific profile and packing primitives; all x86 intrinsics;
- `src/prefix_hash_scan/kernel_common.jl`: geometry-neutral scalar and lookup
  helpers;
- `src/prefix_hash_scan/generic.jl`: compiled motif-specialized scan kernel,
  shared by raw FASTA bytes and decoded 2bit ranges;
- `src/prefix_hash_scan/verification.jl`: Myers, traceback, and result commit;
- `src/prefix_hash_scan/streaming.jl`: chunk streaming and global scheduler;
- `src/prefix_hash_scan/twobit.jl`: 2bit metadata and range decoding.

## Core idea

The algorithm moves part of the alignment work to the query side:

1. Enumerate every symbolic way in which a 16-base reference prefix can be
   produced from a guide while spending at most the requested edits.
2. Apply those symbolic paths to each concrete guide and encode the resulting
   16-mers as 32-bit integers.
3. Scan only geometry-compatible genome windows and test whether their 16-mer is
   in the guide-derived set.
4. Run exact edit-distance verification only for guide/window pairs that pass
   the prefix test.
5. Compute a traceback only for verified off-targets.

The symbolic paths do not contain completed alignments. They describe possible
prefix alignment topologies, including paths caused by substitutions,
insertions, and deletions. Applying a path to a concrete guide produces a
concrete 16-base hash. This is why the method filters edit-distance candidates,
not only Hamming-distance candidates.

The prefix test is necessary but not sufficient. A hash hit is a candidate;
Myers and the final traceback establish the full edit distance.

## Geometry and kernel

`resolve_prefix_scan_geometry` returns a `PrefixScanGeometry{Kind,M}` for every
eligible motif, or `nothing`. It supplies guide, PAM, prefix, distance,
candidate-span, overlap, and compiled motif matching to validation and
orchestration. `Kind` is a label (`:cas9`, `:cas12a`, `:generic`) used only by
the AVX-512 `:auto` policy and statistics. Kernels dispatch on the matcher type
`M`.

`resolve_generic_prefix_scan_geometry` converts motif properties into the
`PrefixScanMatcher` type: enabled strands, constrained IUPAC positions, guide
offsets after PAM removal, orientation, and coordinate offsets. Generated
functions turn that type into straight-line validity, PAM-matching, and prefix
packing operations. The 64-start SIMD block loop contains no runtime motif
branches.

Examples:

- Cas9: 23-base window, `20N + NGG`, extends in the 5-prime direction;
- Cas12a: 25-base window, forward `TTTV + 21N`, reverse `21N + BAAA`, with the
  opposite extension and coordinate rules;
- a 25-base guide followed by an `NNT` PAM resolves to
  `PrefixScanGeometry{:generic}` with a 28-base candidate span:

```julia
motif = Motif(
    "25N_NNT",
    repeat("N", 25) * "NNT",
    repeat("X", 25) * "NNT",
    true, true, 4, true, 0,
)
```

An all-`X` PAM description produces an empty PAM range and therefore a PAMless
search.

## Backends

### Scan backend

`search_prefixHashScan(...; scan_backend=:auto)` is resolved by
`resolve_prefix_hash_scan_plan`:

| Backend | Selected by `:auto` when | Reference access |
|---|---|---|
| `:streaming_fasta_simd` | eligible geometry, FASTA reference | raw FAI range reads |
| `:streaming_2bit_simd` | eligible geometry, `.2bit` reference | decoded 2bit ranges |
| `:legacy` | no eligible geometry, including `hash_len` other than 16, or `query_variant=:bruteforce` | `findguides` over converted `LongDNA` chromosomes |

A pinned streaming backend that does not match the reference format is an
error. Each query holds at most 64 guides in a `UInt64` mask; the public API
splits larger lists into batches, and the low-level four-argument method
rejects them.

### SIMD backend

`simd_backend=:auto` selects AVX-512 only for benchmark-qualified CPU-family and
geometry-label pairs (`:cas9`, `:cas12a`). Other geometries use AVX2, and CPUs
without AVX2/BMI2 use `:portable`. `:avx2`, `:avx512`, and `:portable` force a
backend; `:avx2` and `:avx512` error when the CPU lacks them.

`:portable` builds the same 64-bit `A`/`C`/`G`/`T` profiles from eight
little-endian 64-bit loads. Each byte is case-folded with `0xdf` and compared
with SWAR: the high bit of a byte is set exactly when the byte equals the
pattern, and one multiplication gathers the eight high bits. Prefix packing
spreads the low and high base bits with shifts and masks instead of `PDEP`. The
compiled kernel contains no x86 intrinsics and no target features, which
`test/src/prefix_hash_scan.jl` and `scripts/verify_simd_codegen.jl portable`
check.

### Reference variants

Tests compare the streaming backends against independent references:

- `scan_backend=:legacy` uses `findguides`, a `Dict` query, its own prefix
  extraction (`normalized_candidate_prefix`), and `align`, so it shares no
  code with the generic geometry or kernel;
- `query_variant=:bruteforce` skips the prefix filter and verifies every motif
  candidate on the legacy engine, which checks the filter for false negatives;
- the allocating `scan_generic_prefix_hits_raw_range` wrapper performs the
  genome-order (non-bucketed) directory lookup.

Equivalent optimization-era variants were removed; see
[Appendix B](#appendix-b-removed-experimental-backends).

## Pipeline

### 1. Load symbolic prefix paths

`load_prefix_hash_scan_paths` loads exact precomputed paths by guide length,
distance, and prefix length. Existing 20-base and 21-base assets are reused
regardless of PAM sequence or position. If an asset is unavailable, as for a
25-base guide, paths are generated once before guide batching with
`build_PathTemplates`, restricted to the requested prefix and distance,
deduplicated, and reused by every batch. Path generation does not select the
scan kernel.

At p16, each motif has 1, 129, 7,873, 302,337, and 8,196,801 distinct symbolic
paths for d0 through d4 respectively. Paths are shared by every guide in the
query. The d4 matrices occupy about 125 MiB.

### 2. Build concrete hashes for each guide

Guides are oriented to match the prefixHashDB template convention for the
motif's extension direction and converted to a two-bit alphabet. For every
symbolic path, `fill_prefix_hashes_columnwise!` selects 16 positions from the
formatted guide and folds them into a `UInt32`:

```text
hash = 0
for symbolic_position in path
    hash = (hash << 2) | guide_base[symbolic_position]
end
```

Hashes are sorted and deduplicated per guide. Formatting, folding, sorting, and
deduplication run as bounded tasks, at most `scan_threads` workers, each writing
a distinct guide-list slot. `:auto` uses serial construction for one guide or
one worker.

### 3. Merge hashes into the compact query directory

One bit in a `UInt64` identifies each of at most 64 guides. Equal hashes from
different guides are merged by a serial, deterministic heap merge:

```text
concrete 16-mer hash -> 64-bit mask of compatible guides
```

The lookup is not a Julia `Dict`. Sorted 32-bit hashes are split into:

- a direct bucket-offset array selected by the high hash bits;
- compact `UInt16` suffixes inside each bucket;
- a parallel `UInt64` guide-mask array.

A 26-bit presence bitmap is checked first. Most genome hashes fail this test and
never access the larger directory. The bitmap may create extra work, but cannot
reject a hash that exists in the directory.

### 4. Stream reference chunks

Work is scheduled globally as `(chromosome, 2 MiB core range)` items claimed
through an atomic counter. Each worker keeps its own reference handle, byte
buffer, plus/minus `PrefixHashScanHit` pair, and candidate/radix scratch
buffers. Each read includes a left edit-distance overlap and a right
candidate-span/extension overlap from the geometry; emission is restricted to
the core range.

FASTA chunks are read through the `.fai` offsets and line geometry, and
newlines are removed in place. 2bit chunks are decoded from packed bases and
N-block metadata into the same ASCII buffer. The entire chromosome is never
converted to `LongDNA`.

### 5. Detect motif windows with SIMD

`scan_generic_prefix_hits_raw_range!` (and its `_bucketed!` variant) first
clears the worker hit vectors, then profiles raw ASCII reference bytes in
blocks. AVX2 uses two 32-byte loads, AVX-512BW one 64-byte load with mask
comparisons, and `:portable` eight 64-bit SWAR loads. All produce identical
64-bit `A`, `C`, `G`, and `T` profiles. Adjacent profiles form a 128-base view,
enough to evaluate 64 candidate starts together.

The generated matcher computes, with bit operations:

- starts whose complete window has at most `ambig_max` ambiguous bases;
- forward-strand starts whose PAM matches;
- reverse-strand starts whose reverse-complement PAM matches.

Only set bits are visited. The oriented 16-base prefix is packed into a
`UInt32` with BMI2 `PDEP` on x86 backends, or shifts and masks on `:portable`.
Reverse-strand hashes are reversed and complemented with bit operations. When a
strand's 16 prefix offsets form one descending run, the two reversals cancel and
the prefix packs directly.

### 6. Apply the symbolic-prefix filter

Presence-bitmap checks run in genome order. Survivors are packed as a `UInt64`
of the 32-bit hash and the local candidate start. Three stable worker-scratch
radix passes order the hash by its 10-bit suffix and two 11-bit bucket digits.
Directory lookup then walks hashes in bucket order and reuses the mask for
repeated hashes.

Accepted hits are sorted back by candidate start before verification, preserving
chromosome, strand, and coordinate order. A hit produces:

```text
(candidate start, mask of potentially matching guides)
```

### 7. Verify full edit distance with Myers

Every candidate bit is verified with `prefix_hash_scan_raw_myers_distance`
against a per-guide equality profile from
`build_prefix_hash_scan_myers_profiles`. The verifier reads the reusable raw
buffer directly, handles both strands, evaluates the guide against the extended
reference, and returns a value above the requested distance for rejected
candidates. It is allocation-free and a full Levenshtein filter, including
indels. Ambiguous reference bases outside the prefix use IUPAC-aware masks.

### 8. Trace back accepted candidates

Only candidates within the threshold are materialized as `LongDNA`. `align`
then computes the alignment strings and distance used in the output, preserving
prefixHashDB-compatible reporting. Count mode skips this step and uses the raw
Myers distance.

### 9. Commit results deterministically

Results are stored at stable work indices and committed in reference order: for
each chromosome, all plus-strand chunks, then all minus-strand chunks. The main
task deduplicates complete output records, applies early-stopping counters, and
writes CSV rows.

### 10. Early stopping

Finite limits activate chunk-local guide counters. A guide retires only when a
`(limit + 1)`th accepted hit proves one exact-distance bucket incomplete.
Workers mask retired guide bits before Myers verification and stop claiming
chunks when all guides are retired. Count output caps each bucket and reports
`complete=false`; non-triggering buckets may then be partial lower bounds.
Detail output keeps any valid capped subset, which may vary with scheduling.

### 11. Progress and memory reporting

`verbose=true` logs these `@info` records. The CLI always sets it.

| Record | When | Fields |
|---|---|---|
| `prefixHashScan guide batching` | once, more than 64 guides | guides, batch size, batches |
| `prefixHashScan execution` | once | geometry, reference format, scan and SIMD backends, scheduler, threads, chunk size |
| `prefixHashScan memory` | each batch, after query construction | `paths`, `transient_guide_hashes`, `query`, `query_build_s` |
| `prefixHashScan progress` | at most every `progress_interval` seconds (default 30) during the scan | `batch`, `percent`, `elapsed_s`, `eta_s` |
| `prefixHashScan batch done` | each batch | `elapsed_s`, `scan_s`, `peak_rss`, `chunks_claimed` |

Memory fields:

- `paths`: the symbolic path matrix. Every batch shares it.
- `transient_guide_hashes`: the sorted per-guide hash lists before the
  cross-guide merge, 4 bytes per hash. They are freed after the merge, so this
  is the main transient peak at d4 (about 450 MB for 61 Cas9 guides).
- `query`: the compact directory and the 8 MiB presence bitmap. The `:legacy`
  engine reports its `Dict` size and no path or guide-hash memory.
- `peak_rss`: `Sys.maxrss()`, the peak resident set of the whole process since
  it started, not of this search.

Progress counts reference bases in finished chunks (streaming) or chromosomes
(`:legacy`). Workers add finished bases to an atomic counter. The worker that
wins a compare-and-swap on the next report time writes the record, so reporting
needs no polling task and adds no latency at the end of the scan. Searches
shorter than the interval log no progress record. The CLI sets it with
`--progress_interval`. `scan_s` includes result commit.
`chunks_claimed` is `"all"` without early stopping and `claimed/total` with it.

`PrefixHashScanStats` stores the same sizes in `path_bytes`,
`guide_hash_bytes`, `query_bytes`, and `peak_rss_bytes`. With `verbose=false`
the only added work is one `nothing` check per chunk or chromosome.

## Pseudocode

```text
geometry = resolve_geometry(motif, distance=k, prefix=16)
paths = load_or_generate_symbolic_paths(geometry, distance=k, prefix=16)

for batch in partition(guides, 64):
    for guide in batch (parallel):
        hashes[guide] = unique(sort(apply_each_path(paths, orient(guide))))
    query = compact_directory(heap_merge_into_guide_masks(hashes))
    query = add_presence_bitmap(query, bits=26)
    myers_profiles = build_myers_profiles(batch)

    parallel workers claim (chromosome, 2 MiB chunk) items:
        bases = read_chunk_with_overlap(chunk)        # FASTA or 2bit
        pam_masks = generated_simd_match(matcher_type, bases)

        for candidate_start in set_bits(pam_masks):
            hash = pack_oriented_16mer(bases, candidate_start)
            if presence_bitmap_contains(query, hash):
                append_packed_candidate_by_strand(hash, candidate_start)

        for strand_candidates in plus then minus:
            for hash, candidate_start in radix_order(strand_candidates):
                guide_mask = directory_lookup_or_reuse(query, hash)
                if guide_mask != 0:
                    append_hit(strand_hits, candidate_start, guide_mask)
            sort_hits_by_candidate_start(strand_hits)

            for hit in strand_hits:
                for guide in set_bits(hit.guide_mask & active_guides):
                    d = raw_myers_distance(guide, bases, hit.candidate_start)
                    if d <= k:
                        detail: retain(traceback(guide, materialize(hit)))
                        counts: increment(guide, d)

    commit_in_reference_order()
```

## Comparisons

### Versus the general paths

| Stage | Streaming backends | `:legacy` |
|---|---|---|
| Query structure | Presence bitmap + compact directory + `UInt64` guide masks | `Dict` from hash to `UInt64` guide mask |
| Reference access | Raw FASTA/2bit chunk reads into reusable buffers | Converted `LongDNA` chromosomes |
| Motif search | Generated matcher on 64-start SIMD blocks | `findguides` |
| Distance rejection | Allocation-free raw Myers before materialization | `align` |
| Parallelism | Global chunk scheduler | One record at a time |
| Purpose | Production | Configurations outside the envelope; correctness reference |

### Why it can beat prefixHashDB

`prefixHashDB` does not scan the genome at search time, but it must load and
traverse a large genome-derived index with many partitions and locations.
`prefixHashScan` instead performs a predictable sequential reference pass and
keeps the smaller guide-derived query structure in memory. On this machine and
workload, sequential scanning plus SIMD filtering is cheaper than the
prefixHashDB search-time index access.

This does not mean an index is intrinsically slower. Results depend on storage,
cache state, guide count, index layout, and whether index construction and
storage are included. The architectural cost model is:

```text
prefixHashDB total = database build + searches * indexed search
prefixHashScan total = searches * ceil(guides / 64) * reference scan
```

Thus 61, 1,024, and 4,096 guides require 1, 16, and 64 reference scans. A
persistent index can win for repeated searches even when one scan is faster than
one indexed search.

### Why it differs from Sassy

Sassy performs SIMD Myers-style approximate matching while scanning the text.
`prefixHashScan` uses SIMD mainly to classify bases, identify PAM windows, and
pack exact 16-mer equality keys. Symbolic prefix paths reject most guide/window
pairs before full Myers verification.

The tradeoff is memory:

- Sassy stores compact pattern state and calculates more alignment state during
  the scan.
- `prefixHashScan` materializes millions of guide-derived prefix hashes and
  performs random presence checks followed by bucket-ordered directory lookups.

Which approach wins depends on query-table size, cache behavior, number of
guides, PAM density, and candidate rate. On the PAMless GRCh38 workload, where
Sassy needs no PAM post-filter, prefixHashScan was 5.4x to 7.3x faster than
the fastest Sassy configuration at d0 through d3, and 1.4x faster at d4
([PAMless benchmark](#pamless-versus-cas9-and-sassy)).

## Current performance

GRCh38, 61 guides per motif. These numbers are evidence for one shared host and
workload, not a universal speed claim. Setups and older records are in
[Appendix A](#appendix-a-benchmark-history).

Distance sweep, detail output, 24 threads, July 22, 2026 (prefixHashDB build
excluded):

| Motif | Distance | Results | `prefixHashDB` median | `prefixHashScan` median | Scan speedup |
|---|---:|---:|---:|---:|---:|
| Cas9 | 0 | 65 | 14.661 s | 0.618 s | 23.7x |
| Cas9 | 1 | 130 | 14.732 s | 0.686 s | 21.5x |
| Cas9 | 2 | 1,244 | 15.251 s | 0.742 s | 20.6x |
| Cas9 | 3 | 25,826 | 22.355 s | 1.875 s | 11.9x |
| Cas9 | 4 | 381,003 | 136.837 s | 13.131 s | 10.4x |
| Cas12a | 0 | 8,173 | 7.326 s | 0.499 s | 14.7x |
| Cas12a | 1 | 26,612 | 8.309 s | 0.574 s | 14.5x |
| Cas12a | 2 | 95,531 | 9.268 s | 1.229 s | 7.5x |
| Cas12a | 3 | 364,581 | 14.241 s | 3.532 s | 4.0x |
| Cas12a | 4 | 1,073,287 | 120.346 s | 18.451 s | 6.5x |

The d4 prefixHashDB builds took 1,419.7 s for Cas9 and 634.6 s for Cas12a and
occupied 4.60 GB and 2.95 GB.

Count versus detail output, 24 threads, August 8, 2026:

| Motif | Distance | Detail median | Count median | Speedup |
|---|---:|---:|---:|---:|
| Cas9 | 3 | 1.556 s | 1.446 s | 1.08x |
| Cas9 | 4 | 14.376 s | 12.556 s | 1.14x |
| Cas12a | 3 | 3.819 s | 1.069 s | 3.57x |
| Cas12a | 4 | 20.819 s | 10.758 s | 1.94x |

SIMD backends, October 5, 2026, Xeon Gold 6126 (Skylake-SP):

| Measurement | AVX2 | AVX-512 | Portable | `:fused_directory` |
|---|---:|---:|---:|---:|
| Cas9 scanner, 32 MB, 1 thread | 68.0 ms | 65.0 ms | 68.6 ms | - |
| Cas12a scanner, 32 MB, 1 thread | 29.0 ms | 25.5 ms | 30.2 ms | - |
| Cas9 GRCh38 d3 detail, 24 threads | 1.359 s | - | 1.365 s | 11.481 s |
| Cas12a GRCh38 d3 detail, 24 threads | 3.557 s | - | 3.959 s | 14.647 s |

End-to-end values are medians of 5 alternating runs; outputs were identical.
`:fused_directory` was removed after this measurement (Appendix B).
On x86 the portable profile costs little because LLVM vectorizes the SWAR loop
and the scan is a small share of the search. ARM performance is not measured.

### PAMless versus Cas9 and Sassy

GRCh38, the same 61 Cas9 guides (20 nt, no PAM), 24 threads, unlimited early
stopping, October 6, 2026. PAMless uses `pamless=true`; Cas9 uses `NGG`. Sassy
is stock `sassy search` v0.2.6 (master `7c9f8fc`) with `--max-n-frac 0` to
match `ambig_max=0`. Sassy always computes tracebacks, so its closest CHOPOFF
equivalent is detail output; d3 and d4 use count output because PAMless
detail would be multi-GB. Medians of 5 runs (d0-d3, prefixHashScan), 3 runs
(d4, prefixHashScan), 3-4 runs (Sassy d0-d3), and 1 run (Sassy d4). The host
was shared and loaded (load average 16-60 on 48 cores), so expect
run-to-run variation of about 30%.

| Distance | Output | PAMless | Cas9 | PAMless/Cas9 | Sassy `-a dna` | Sassy IUPAC (default) | Sassy/PAMless |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | detail | 3.98 s | 0.64 s | 6.2x | 22.6 s | 60.7 s | 5.7x |
| 1 | detail | 4.95 s | 0.61 s | 8.1x | 26.5 s | 41.0 s | 5.4x |
| 2 | detail | 4.67 s | 0.69 s | 6.7x | 34.2 s | 48.1 s | 7.3x |
| 3 | counts | 6.37 s | 1.09 s | 5.8x | 39.7 s | 52.8 s | 6.2x |
| 4 | counts | 48.97 s | 15.42 s | 3.2x | panics | 67.0 s | 1.4x |

The Sassy/PAMless ratio uses the faster Sassy configuration that completed.
Result rows (PAMless detail; Sassy): d0 75; 75. d1 848; 666. d2 23,792;
20,236. d3 Sassy 408,883 and d4 Sassy 5,182,770.

Counters: PAMless scans 5,891,652,878 candidate windows, 19.4x the
304,418,266 Cas9 windows. Prefix hits grow by 17.8x at d2 (1,571,505 versus
88,450), 16.2x at d3, and 19.3x at d4. The cost of losing the PAM filter is
6-8x at d0-d2, not 19x, because the query and its construction are the same
size and the scan kernel is not the only cost. d4 is dominated by query
construction (12.3 s PAMless, 7.4 s Cas9, measured in a stats pass) and by
280 million guide/window verifications.

Correctness:

- Every Cas9 detail row at d0-d2 appears in the PAMless output with the same
  alignment and distance. `start` shifts by -3 on `+` and +3 on `-` (1,439 of
  1,439 rows).
- Every Sassy match at d1 and d2 has a PAMless row with the same distance
  within ±k bases (666/666 and 20,236/20,236). At d0, per-guide counts are
  identical.
- PAMless rows exceed Sassy matches because CHOPOFF reports end-gap shadow
  alignments, for example `...GGGG-C` one base next to an exact hit, as
  separate loci. Sassy reports one local minimum per locus. Every PAMless row
  at d1 and d2 lies within ±k of a Sassy match of equal or lower cost, except
  4 (d1) and 42 (d2) shadows 2 bases from a Sassy match.

Sassy observations:

- Sassy parallelizes across FASTA records, not within a record, so chr1 sets
  its wall time. prefixHashScan schedules 2 MiB chunks globally.
- With the default IUPAC alphabet, `N` matches every base. N blocks generate
  hits that are filtered only afterwards: one all-N 5 Mb chr21 slice took
  6.4 s with v1 and 60 s with v2. `-a dna` avoids this (0.6 s) but panics at
  d4 when a traceback reaches an `N` (`src/trace.rs:376`), and gives one extra
  row at d3.
- v2 (`search --v2`) is about 4x faster than v1 on clean sequence, but it took
  582 s on GRCh38 at d2 because of N blocks. `--v2 -a dna` returned no matches.
  The `crispr_v2_search` branch (March 2026, 133 commits behind master) adds
  only `sassy crispr --v2`. v2 search is already merged into master.

Reproduce with `test/local_human/benchmark_human_pamless.jl` and
`test/local_human/run_sassy.sh`.

Workload shape, Cas9 versus Cas12a at d3:

| Counter | Cas9 | Cas12a |
|---|---:|---:|
| GRCh38 motif candidates | 304,418,266 | 136,277,556 |
| Exact directory hits | 1,547,796 | 1,546,515 |
| Guide/window verification pairs | 1,583,279 | 1,887,411 |
| Tracebacks and emitted rows | 25,826 | 364,581 |
| Precomputed path rows | 302,337 | 302,337 |
| Concrete query hashes before cross-guide merge | 7,044,938 | 6,947,869 |

Cas9 time is dominated by the scan/lookup kernel. Cas12a performs only 19% more
verifications but 14.1x more tracebacks and output commits, so its detail time
is dominated by materialization, `align`, deduplication, and CSV commit.

D4 is a query-construction problem as well as a scan problem: Cas9 d4 expanded
8,196,801 paths into 111,720,240 guide/hash associations and a 1.09 GB query,
with 7.61 s query construction.

## Limitations

1. Each optimized query holds at most 64 guides. Larger lists rescan the
   reference once per batch.
2. The prefix length is fixed at 16. Other prefix lengths and motifs outside
   the envelope use the slower `:legacy` engine.
3. Query guides must be unambiguous. `ambig_max` above 3 is unsupported.
4. ARM performance of `:portable` is unmeasured; it is verified only by its
   target-independent IR.
5. Early stopping cannot cancel chunks already claimed by workers. Chunk-local
   reduction bounds this overshoot.
6. Global scheduling uses concurrent reference seeks. Evidence is for
   warm-cache GRCh38; cold-cache and networked filesystems are unmeasured.
7. Buffered workers retain the largest observed per-chunk hit and
   prefilter-survivor capacities, which raises cumulative allocation by
   1.8-3.4%.
8. Requested statistics add counters and `time_ns()` calls to hot loops;
   `stats=nothing` compiles them out.
9. `PrefixHashScanStats` mixes summed worker CPU times with wall-clock fields;
   they cannot be compared directly.
10. Query construction is rebuilt for every call. Per-guide hash lists are
    parallel, but the heap merge and directory construction are serial.
    Cross-run reuse is out of scope for the one-shot workload.
11. The 8.4 MB presence bitmap is probed in genome order and may be
    memory-latency bound.
12. Progress counts only finished work. When early stopping retires every
    guide, the last progress record stays below 100%. `peak_rss` is a
    process-lifetime peak, so a later batch cannot report a lower value than an
    earlier one.
13. Guide lengths above 28 bases are not qualified at d4: prefixHashDB, the
    current oracle, uses a packed representation limited to 32 bases.

## Roadmap

`prefixHashScan` is intended to become the primary CHOPOFF algorithm for
ordinary reference-genome search: exact, deterministic,
prefixHashDB-compatible coordinates and distances, no genome-specific database,
and one bounded reference scan per 64-guide batch.

### Completed

- Registered and custom motifs at d0 through d4 in the Julia API and the
  standalone CLI, including PAM-left, PAM-right, internal-PAM, PAMless,
  strand-subset, and extension-direction definitions. `verbose` reports the
  resolved backend without enabling statistics.
- One generic motif-specialized kernel for all eligible motifs; the hand-written
  Cas9 and Cas12a kernels were removed on October 4, 2026.
- Distance 4 at p16 as the functional ceiling. The p14/p15/p16 evaluation
  selected p16.
- Sequential 64-guide batching for larger lists, with atomic multi-batch detail
  output.
- Count output without traceback, with exact parity against
  `summarize_offtargets(detail; distance=k)` as a release gate.
- Computational early stopping that masks retired guides and cancels future
  chunk claims.
- 2bit streaming and bounded IUPAC reference ambiguity (`ambig_max=0:3`).
- AVX-512F/BW and `:portable` backends with parity tests and codegen
  verification.
- Progress and path, guide-hash, query, and peak-RSS memory reporting through
  `verbose` and `PrefixHashScanStats` (October 6, 2026).
- Representative generic qualification against prefixHashDB: full GRCh38
  Cas9-NGA, CasX, and 25-base-guide cases, 65-guide multi-batch searches, d0
  through d4, ambiguity zero through three, and bounded internal-PAM, PAMless,
  16-base-guide, strand-subset, FASTA/2bit, IUPAC, indel, and chunk-boundary
  cases. All 220 detail cases passed the reference-backed parity classifier and
  all 220 count comparisons passed. Of the detail cases, 147 were exact; the
  remainder contained only classified prefixHashDB ambiguity-limit or
  duplicate-row behavior.

### Closed decisions

- Distance 4 stays at p16 and is a stretch configuration, not an optimization
  target. Compressed or staged d4 representations are out of scope.
- Sequential 64-guide batching is the large-guide architecture. Reconsider a
  one-pass design only if a bounded-memory version beats batching by at least
  10% on both dispersed and related guide sets.
- No automatic crossover policy between prefixHashScan and prefixHashDB is
  planned.

### Remaining work, in priority order

1. **Qualification maintenance.** Add randomized property tests and qualify
   guide lengths above 28 with an oracle other than prefixHashDB.
2. **Portable performance.** Measure `:portable` on ARM (for example Graviton or
   Apple Silicon) and on AMD Zen1/Zen2, where microcoded `PDEP` may make
   portable packing faster than `:avx2`. Add an ARM SIMD path only if profiling
   justifies it.

### Replacement gates

Make `prefixHashScan` the documented default for one-shot searches with eligible
motifs when:

1. Canonical searches retain exact detail parity, and representative generic d0
   through d4 retain exact or reference-classified parity on sample and
   human-scale fixtures.
2. Legacy, FASTA streaming, and 2bit streaming backends have identical
   results under every SIMD backend;
   unsupported configurations select a correct fallback rather than fail.
3. Count output marks early-stopped rows incomplete; detail output documents its
   scheduling-dependent valid subset.
4. Progress/memory reporting is clear enough for production use.
5. The `:portable` backend is qualified on ARM hardware.
6. Cas9/d3 and Cas12a/d3 detail latency regresses by no more than 3% unless a
   measured feature-level benefit justifies it.

Gates 1-4 and 6 are met. Keep `prefixHashDB` as the persistent-index backend
for repeated or heavily capped workloads.

## Performance research

Speed research is optional and does not precede the remaining work above.
Continue an experiment only with exact parity and a credible route to at least a
10% end-to-end improvement on GRCh38 in its intended workload. Report detail and
count modes separately, and verify that a motif-specific win does not regress
the other motif.

Another 1.5-3x is plausible only if profiling confirms avoidable lookup,
scheduling, or temporary-data costs. Another 10x is unlikely because every exact
indexless search must still inspect the reference.

Open directions:

- **Cas9:** genome-order presence-bitmap latency and last-level-cache behavior;
  compare the bitmap/directory with a blocked Bloom, xor, or quotient-style
  prefilter at equal memory; NUMA placement, pinning, and per-socket query
  replicas; scheduler tail imbalance.
- **Cas12a detail:** SIMD-vectorized Myers or traceback across independent
  candidates; avoiding string conversion and generic dedup keys during commit.
- **I/O:** fewer FASTA newline-compaction copies or a direct mmap-backed scan,
  only after I/O measurements.

Open questions:

- How much of the 0.32-0.34 s query build is the serial heap merge and directory
  construction?
- Are the genome-order presence-bitmap probes limited by last-level cache misses
  or memory latency?
- How much of the 24-core 44.4% wait share is tail imbalance versus task/runtime
  overhead, and can NUMA-local query placement reduce it?
- Why does immediate verification show no latency gain despite allocating 5.78%
  fewer bytes: instruction pressure or phase-locality loss?
- How much Cas12a time can be removed by batching traceback or replacing
  string-based commit/deduplication while preserving exact output ordering?
- Why are Cas12a query/path preparation costs higher than Cas9 despite nearly
  identical concrete-hash and path counts?
- Does Cas12a keep its 2.9x advantage for guide sets with fewer accepted off-targets,
  where traceback and CSV output do not dominate?

## Appendix A: Benchmark history

Entries are newest first. Unless stated otherwise: GRCh38, 61 guides per motif,
warm cache, exact output parity in every comparison. The host was shared, and
identical runs varied by up to 40%.

### PAMless and Sassy (October 6, 2026)

Results are in [PAMless versus Cas9 and Sassy](#pamless-versus-cas9-and-sassy).
A first 3-run pass at load average about 55 gave PAMless d0-d2 medians of
4.7-7.8 s with no trend by distance. A 5-run rerun at load 16-33 gave the
reported d0-d3 values. A chr21 pilot (61 guides, k=3) measured Sassy v1 at
15.2 s and v2 at 65.9 s, which exposed the N-block cost. Sassy was built with
`RUSTFLAGS="-C target-cpu=native"` and Rust 1.99.0. The July
`run_rust_sassy_v2.sh` baseline was removed. It ran a custom
`chopoff_batch_crispr` binary with CHOPOFF-side PAM filtering, which is no
longer in the Sassy source tree. `run_sassy.sh` replaces it with stock
`sassy search` and `sassy crispr`. Outputs:
`test/local_human/outputs/pamless_20261006*`.

### Generic kernel replaces hand-written kernels (August 2 to October 4, 2026)

An August 2 scanner microbenchmark measured the generic kernel at 1.178x the
latency of the hand-written Cas9 kernel. Two generator changes closed the gap
before the hand-written kernels were removed:

- strands whose 16 prefix offsets form one descending run no longer emit
  `bitreverse` on both profiles followed by `reverse_codes`;
- Myers and materialization reference offsets fold to
  `base + step * ref_idx` when the guide and extension form one run, instead of
  a tuple lookup per base.

Generic latency relative to hand-written code, one thread, 32 MB seeded random
A/C/G/T, d3, 11 alternating runs after 2 warmups:

| Stage | Cas9 before | Cas9 after | Cas12a before | Cas12a after |
|---|---:|---:|---:|---:|
| Raw scan, AVX2 | 1.164x | 0.995x | 1.114x | 0.995x |
| Raw scan, AVX-512 | 1.295x | 0.961x | 1.126x | 0.919x |
| Raw Myers verification | 1.042x | 1.005x | 1.103x | 0.986x |
| Raw materialization | 1.6x | 0.81x | 1.5x | 0.96x |
| LongDNA scan (fused/legacy) | 2.6x | 0.08x | 2.0x | 0.05x |
| LongDNA materialization | 6.6x | 0.71x | 8.6x | 0.79x |

The LongDNA scan now builds profiles from packed 4-bit words and runs the raw
block kernel. End-to-end d3 detail at 24 threads, two rounds of 5 runs: Cas9
1.359 / 1.463 s before versus 1.338 / 1.471 s after; Cas12a 3.793 / 4.146 s
versus 4.435 / 3.712 s. The differences are within host noise. Output was
identical (25,826 and 364,581 rows).

### AVX-512 end-to-end qualification (September 26-27, 2026)

`scripts/benchmark_prefix_hash_scan_avx512.jl`, d3, count output, 11
alternating runs after 2 warmups. AVX-512 end-to-end speedup over AVX2:

| Code | Threads | Cas9 | Cas12a | Cas9_NGA (not `:auto`) |
|---|---:|---:|---:|---:|
| `846c0c17` (before refactor) | 24 | 0.93x | 0.96x | 0.97x |
| refactored | 24 | 0.96x | 0.85x | 0.96x |
| refactored | 12 | 0.99x | 0.99x | 1.02x |
| refactored | 8 | 0.94x | 1.02x | 0.98x |

The gate requires at least 0.97x for `:auto`-eligible motifs. The scanner-only
stage was 1.08x to 1.49x faster with AVX-512 in every run, but the scan is only
about 0.1 s of about 1.5 s end to end, so a 20-30% faster kernel changes
end-to-end time by about 2%, below host noise. The pre-refactor commit fails the
same gate, and there is no thread-count trend. Decision: `simd_backend=:auto`
is unchanged; repeat on an idle host before changing it.

### Count output (August 8, 2026)

24 threads, one warmup, three alternating timed runs per mode. Count/detail
summaries had exact parity and count mode performed zero tracebacks. Results
are in [Current performance](#current-performance).

### Full distance 0-4 sweep (July 22, 2026)

24 threads, one warmup, five timed repetitions per algorithm and distance,
alternating order, unlimited early-stopping thresholds. One d4 prefixHashDB per
motif was built with one thread and reused for all distances. A separate untimed
pass collected phase counters. Results are in
[Current performance](#current-performance).

Raw prefixHashDB detail output is the gold standard. Parity compares exact
detail-row multisets, including duplicate multiplicity. The July 23 parity
repair reports PASS with zero scan-only and zero prefix-only rows for all ten
motif/distance cases. An initial Cas12a `parity=false` came from an obsolete
Cas9-only filter that rejected every `extends5=false` row.

Phase counters: d4 query construction took 7.61 s for Cas9 and 7.26 s for
Cas12a, producing 111,720,240 and 109,038,602 guide/hash associations. Cas12a
performed 364,581 tracebacks at d3 and 1,073,287 at d4. An earlier isolated
61-guide Cas9 d4 run with 24 query workers took 10.8 seconds for query
construction and reached 4.51 GB peak process RSS.

### Matched Cas9 and Cas12a benchmark (July 16, 2026)

d3, 8 threads, unlimited early stopping, existing prefixHashDB indexes. The
Cas12a guides were sampled from 61 distributed canonical `TTTV` sites so every
query had a real on-target. One warmup, three measured runs.

| Motif | Results | `prefixHashScan` median | `prefixHashDB` median | Scan speedup | Scan runs | DB runs |
|---|---:|---:|---:|---:|---|---|
| Cas9 | 25,826 | 2.694 s | 22.861 s | 8.49x | 2.826, 2.694, 2.673 s | 23.263, 22.818, 22.861 s |
| Cas12a | 364,581 | 5.083 s | 14.856 s | 2.92x | 5.377, 5.083, 4.863 s | 15.114, 14.856, 14.583 s |

A single-pass harness measured Cas9 at 4.507 s versus 46.863 s (10.40x) and
Cas12a at 7.921 s versus 19.823 s (2.50x), with more loading noise. The
prefixHashDB indexes occupy about 3.3 GiB (Cas9) and 1.5 GiB (Cas12a). Both
motifs had exact core-result parity on guide, distance, chromosome, position,
and strand.

Sampled profiles: Cas12a adds a large contribution from
`evaluate_prefix_hash_scan_hits!`, `align`, sequence/string materialization,
deduplication, and CSV commit. Thread utilization was about 69% for Cas9 and
71% for Cas12a. Profile `align_ns`, `verify_ns`, and other worker fields are
summed across threads and can exceed wall time. Sampling raised wall time to
5.816 s and 8.853 s, so profiled timings are not used for speedups.

### Large-guide batching (July 2026)

1,024 Cas9/d3 guides, eight threads, seven rotated timed repetitions, exact
detail-row multiset comparison. Brackets are bootstrap 95% intervals for the
median.

| workload | sequential 64-guide batches | best large directory | ratio |
|---|---:|---:|---:|
| dispersed guides | 608 s [556, 628] | 1,810 s [1,722, 1,977] | 2.98× |
| related guides | 407 s [369, 420] | 1,130 s [933, 1,299] | 2.78× |

The rejected one-pass design used a hash-to-guide-ID directory and optional
reference-aware query filtering. Filtering cut the dispersed query from roughly
1.45 GB to 242 MB, but latency stayed about three times slower: tens of millions
of detail rows made one monolithic result and dedup state the dominant cost.
Wider presence filters and a 12-base bucket did not close the gap.

### Early stopping

d3, 24 threads, Cas9 prefixHash-style caps: verification pairs fell from
1,583,279 to 122,692, 51 guides retired, and detail median improved from
1.418 s to 1.359 s. An 11-run single-thread Cas9 d4 count test measured a
2.27x paired speedup (84.94 s versus 36.80 s median). Unlimited and default
one-million caps stayed near baseline.

### Cas9 scaling and tuning record

Warm-cache d3 search on pinned physical CPUs:

| Configuration | 12 CPUs | 24 CPUs | Notes |
|---|---:|---:|---|
| Current: bucketed lookup, parallel query | 1.710 s | 1.067 s | Paired lookup A/B; no statistics |
| Inline-lookup reference, parallel query | 1.833 s | 1.139 s | Same A/B and exact output |
| Previous production confirmation | 1.725 s | 1.200 s | Before bucketed lookup |
| 2 MiB chunks, serial-query reference | 2.064 s | 1.591 s | Historical query A/B |
| 8 MiB chunks, serial-query reference | 2.313 s | 1.717 s | Chunk-confirmation run |
| `prefixHashDB` | 27.237 s | Not measured | Existing index; build excluded |

Each promoted change used 15 alternating pairs unless noted. Every run kept the
same 25,826 results, output bytes, and semantic counters.

| Change | Result | Decision |
|---|---|---|
| Radix-ordered directory lookup | Prepared scan 1.321 → 1.171 s (11.3%) at 12 CPUs, 0.717 → 0.649 s (9.5%) at 24. End-to-end 6.7% and 6.2%. Allocation +1.8% / +3.4%; peak RSS 1.092 → 1.063 GB (2.6%) and 1.362 → 1.259 GB (7.6%). | Promoted |
| Parallel per-guide query construction | Query build 0.728 → 0.320 s (56.1%) at 12 CPUs, 0.778 → 0.343 s (56.0%) at 24. End-to-end 16.4% and 24.6%. 8-guide crossover 5.6% and 8.7%. | Promoted for multi-guide queries |
| 2 MiB chunks (versus 8 MiB) | End-to-end 2.313 → 2.227 s (3.73%) at 12 CPUs, 1.717 → 1.682 s (2.00%) at 24. Allocation 1.010 → 0.869 GB (14.0%) and 1.153 → 0.895 GB (22.3%). Peak RSS 1.650 → 1.393 GB (15.6%). | Promoted |
| Parameter sweep: 2, 4, 8, 16 MiB chunks; 22, 24, 26 prefilter bits; 9 through 12 bucket bases | All configurations preserved output | Kept 26 bits and 11 bases |
| Global chunk scheduling (8 MiB) | Prepared scan 1.593 → 1.407 s (13.21%) at 12 cores; 1.417 → 0.758 s (46.54% lower latency) at 24 cores. 6.24% and 12.62% fewer allocated bytes. | Promoted; whole-chromosome scheduler removed |
| Reused buffered hit vectors | Prepared scan 1.686 → 1.675 s (0.69%); allocated bytes 812,254,224 → 750,851,360 (7.56%) | Promoted |
| Fused scan/verify | 1.610 s fused versus 1.616 s buffered (0.37%, inside variance); allocated bytes 1,113,346,392 → 1,049,027,928 (5.78%) | Rejected, removed |
| No-statistics hot path | 1.09% faster | Promoted |

Reproduce with `scripts/benchmark_prefix_hash_scan_tuning.jl` at commit
`f27312d2` and `CHOPOFF_TUNING_STAGE=chunk|prefilter|bucket|final|query|lookup`,
for example:

```bash
git checkout f27312d2
CHOPOFF_TUNING_STAGE=lookup \
  julia --project=. scripts/benchmark_prefix_hash_scan_tuning.jl
```

The current script keeps only the `chunk` and `final` stages, because the
prefilter width, bucket width, query construction, and lookup order are now
fixed.

Query and profile state at that point: 6,898,183 exact hashes; offset,
suffix/mask, and presence arrays of about 16.8 MB, 69.0 MB, and 8.4 MB. The
26-bit bitmap had 3,121,043 set prefixes (4.65%), implying roughly 14.2 million
prefilter survivors from the 304.4 million PAM candidates. 1,547,796 windows
passed the directory and expanded to 1,583,279 guide/window pairs, 1.023 per
hit.

Production sampling (2 MiB, parallel query, before radix lookup):

| Share of cumulative samples | 12 CPUs | 24 CPUs |
|---|---:|---:|
| SIMD scan/lookup | 35.8% | 16.0% |
| Directory lookup | 19.5% | 7.9% |
| Worker wait | 24.2% | 44.4% |
| Myers verification | 2.6% | 1.8% |
| FASTA range reading | 3.2% | 1.5% |

Hardware performance counters were unavailable. Earlier statistics-enabled
profiling attributed about 0.97 s of a 2.51 s 12-core run to serial query
construction and 1.54 s to scan plus commit, before query construction was
parallelized.

The July 2026 stabilization refactor (source split, `PrefixScanGeometry`)
preserved full-GRCh38 output bytes, signatures, counters, and all 25,826 rows at
12 and 24 CPUs.

## Appendix B: Removed experimental backends

These paths were measured, lost to the current defaults, and were removed. Check
out the listed commit, the last one that contains them, to rerun them. The
measurements in this document are unchanged.

| Removed path | What it did | Why it was removed | Last commit |
|---|---|---|---|
| `scan_backend=:fused_fasta_simd` | Loaded a whole chromosome as raw bytes and split it across threads with the SIMD kernel | Replaced by global 2 MiB chunk streaming: 13.2% faster at 12 cores, 46.5% lower latency at 24 cores, lower peak RSS | `846c0c17` |
| `scan_backend=:streaming_fasta_simd_fused` | Ran Myers verification inside the SIMD scan loop, with no hit vectors | 0.37% latency difference (inside run variance); 5.78% fewer allocated bytes; needed a second copy of each kernel loop and did not support early stopping | `846c0c17` |
| `scan_backend=:fused_dict` | Whole-chromosome fused scan with a Julia `Dict` query | Replaced by the compact bitmap and directory query | `846c0c17` |
| Streaming `Val(:chromosome)` scheduler | One whole chromosome per worker | Tail imbalance; replaced by global chunk scheduling | `846c0c17` |
| `query_variant=:baseline` | Built each guide's hash set one path row at a time | Replaced by `:columnwise`, which produces identical hashes | `846c0c17` |
| Hand-written Cas9 and Cas12a kernels | Motif-specific SIMD scan, Myers, and materialization; removed October 4, 2026 | Replaced by the generic kernel at equal or lower latency | `ef865b99` |
| `scan_backend=:fused_directory` | Whole-chromosome `LongDNA` scan with the compact directory and the generic block kernel | Streaming covers every eligible geometry and was 8.4x faster for Cas9 d3 (1.359 s versus 11.481 s); `:auto` could no longer select it | `f27312d2` |
| `verify_variant` (`:align`, `:distance_first`, `:myers_raw`) | Chose the verifier per backend | Only `:fused_directory` used the choice; streaming always uses raw Myers and `:legacy` always uses `align` | `f27312d2` |
| `lookup_variant=:inline` and stream mode `:buffered_reuse` | Directory lookup in genome order | Radix-ordered lookup improved prepared scans by 9.5-11.3% (Appendix A) | `f27312d2` |
| `prefilter_bits` and `bucket_bases` keywords | Configurable prefilter width (0, 22, 24, 26) and bucket width | The sweep kept 26 bits and 11 bases, and bucketed lookup requires both | `f27312d2` |
| `query_build_backend=:serial` | Built per-guide hash lists on one task | Parallel construction improved query build by 56% with identical output (Appendix A) | `f27312d2` |
| `query_variant=:columnwise` and the `CHOPOFF_PREFIX_HASH_SCAN_QUERY` variable | `Dict` from hash to guide-index vectors for more than 64 guides | Sequential 64-guide batching replaced one-pass large queries (Appendix A) | `f27312d2` |
| `scripts/benchmark_prefix_hash_scan_experiment.jl`, `scripts/profile_prefix_hash_scan_query.jl` | Swept the backend, verify, bucket, prefilter, and query variants above | Nothing left to compare | `f27312d2` |
| `candidate_prefix_hashes_direct` | Legacy-engine prefix hashing through the generic geometry's offsets | Reached only when tests pinned `:legacy` on an eligible motif; it made the legacy reference depend on the code it checks | `f27312d2` |
| `CAS9_D3_PREFIX_SCAN_GEOMETRY`, `CAS12A_D3_PREFIX_SCAN_GEOMETRY`, Cas9-default `stream_prefix_hash_scan` and `prefix_hash_scan_raw_myers_distance` methods, `prefix_hash_scan_guide_hashes`, Dict-input `build_prefix_hash_scan_directory` | Test conveniences inside the package | Tests now define them locally | `f27312d2` |
