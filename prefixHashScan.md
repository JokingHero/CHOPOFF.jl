# prefixHashScan

## Status and scope

`prefixHashScan` is an indexless CRISPR off-target search. Its supported public
entrypoint covers Cas9 and Cas12a at edit distances 0 through 4 with a fixed
16-base prefix; broader tuning remains an experimental benchmark and parity
engine. It reuses the
symbolic prefix paths from `prefixHashDB`, builds a guide-specific query
structure in memory, and scans the reference genome directly.

The current optimized geometries are:

- Cas9: 20-base guide with an `NGG` PAM;
- Cas12a: 21-base guide with a `TTTV` PAM;
- edit distance 0 through 4
- 16-base prefix
- at most 64 guides
- FASTA reference with a standard `.fai` file
- x86 CPU with AVX2 and BMI2

Current production tuning for this path is:

- 2 MiB globally scheduled chunks;
- a 26-bit presence prefilter and 11-base directory bucket;
- parallel per-guide hash construction for multi-guide queries, bounded by
  `scan_threads`;
- the deterministic serial heap merge and buffered scan/verify backend;
- radix-ordered compact-directory lookup for 26-bit prefilter survivors.

`query_build_backend=:serial` remains the query-construction reference.
`lookup_variant=:inline` remains the genome-order lookup reference. One-guide
query construction stays serial. Compatible production searches use bucketed
lookup automatically; multi-guide searches also use parallel query construction.

"Indexless" means that no CHOPOFF genome database must be built. The optimized
FASTA reader still requires the small, standard `.fai` random-access index.

The implementation is split by stable responsibility:

- `src/db_prefix_hash_scan.jl`: shared types, orchestration, and public API;
- `src/prefix_hash_scan/query.jl`: symbolic paths, hashes, directory, prefilter;
- `src/prefix_hash_scan/kernel_common.jl`: geometry-neutral SIMD/lookup primitives;
- `src/prefix_hash_scan/cas9.jl`: scalar and AVX2/BMI2 Cas9 scan kernels;
- `src/prefix_hash_scan/cas12a.jl`: scalar and AVX2/BMI2 Cas12a scan kernels;
- `src/prefix_hash_scan/verification.jl`: Myers, traceback, and result commit;
- `src/prefix_hash_scan/streaming.jl`: indexed FASTA reads and global scheduler.

`PrefixScanGeometry{Kind}` supplies guide, PAM, prefix, distance,
candidate-span, and overlap values to validation and orchestration. Literal
Cas9 and Cas12a constants stay inside separate SIMD hot loops.

## Core idea

The algorithm moves part of the alignment work to the query side:

1. Enumerate every symbolic way in which a 16-base reference prefix can be
   produced from a guide while spending at most the requested edits.
2. Apply those symbolic paths to each concrete guide and encode the resulting
   16-mers as 32-bit integers.
3. Scan only geometry-compatible genome windows and test whether their 16-mer is in
   the guide-derived set.
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

## Optimized geometries

Cas9 uses a 23-base window with an `NGG` PAM. Cas12a uses a separate 25-base
window: forward `TTTV + 21N`, reverse `21N + BAAA`, with opposite extension
and coordinate rules. Both use prefix 16 and shipped precomputed path assets.
Dispatch occurs before the scalar or SIMD hot loop.

### Cas12a specialization

Cas12a is not implemented by substituting constants into the Cas9 SIMD loop.
`resolve_prefix_scan_geometry` selects `PrefixScanGeometry{:cas12a}` for the
canonical 21-base `TTTV` motif at distances 0 through 4 with a 16-base prefix.
The named d3 geometry constant is the historical baseline; resolution constructs
the requested distance-specific geometry before dispatch. Scalar or AVX2/BMI2
functions in `cas12a.jl` execute without a motif-kind branch inside the hot loop.

The Cas12a kernel evaluates 25-base candidate windows. It recognizes forward
`TTTV + 21N` sites and reverse-complement `21N + BAAA` sites, rejects a window
containing non-ACGT reference bases, and packs the correctly oriented 16-base
prefix with BMI2. Its PAM-left geometry extends and normalizes candidates in
the opposite direction from Cas9. Shared orchestration obtains candidate span,
chunk overlap, bounds, and verification orientation from the typed geometry
rather than from Cas9 literals.

Cas12a loads the exact precomputed symbolic paths for the requested distance
(302,337 rows at d3), then uses the compact `UInt64` guide-mask query, 26-bit
presence prefilter, radix/bucket directory, global FASTA chunk scheduler, raw
Myers rejection, accepted-hit traceback, and deterministic commit pipeline
described below. `scan_backend=:auto` selects `:streaming_fasta_simd` under the
same FASTA, guide-count, AVX2, and BMI2 requirements as Cas9. The public API
selects it with `motif="Cas12a"`; the CLI uses `--motif Cas12a`.

## Optimized streaming path

### 1. Select the backend

`search_prefixHashScan(...; scan_backend=:auto)` selects
`:streaming_fasta_simd` when all of the following hold:

- the motif has the canonical Cas9 or Cas12a guide/PAM geometry;
- the distance is between 0 and 4;
- `hash_len == 16`;
- the query uses the 64-bit guide mask representation;
- there are no more than 64 guides;
- the reference is FASTA;
- AVX2 and BMI2 are available.

Canonical Cas9 and Cas12a geometries use the same backend-selection policy. If
raw FASTA SIMD is unavailable, `:auto` selects the intermediate
`:fused_directory` backend. Unsupported geometries select `:legacy`.

`:streaming_fasta_simd_fused` is an explicit experimental backend. It verifies
nonzero guide masks immediately in the SIMD scan loop and avoids
`PrefixHashScanHit` vectors. The buffered `:streaming_fasta_simd` backend remains
the `:auto` choice because fusion showed no measurable GRCh38 latency advantage.
Buffered workers retain and clear one plus/minus hit-vector pair across chunks;
the allocating scanner wrapper remains the parity reference.

Streaming work is scheduled globally as `(chromosome, 2 MiB core range)` items.
Workers claim items atomically and keep independent FASTA handles, read buffers,
and hit scratch vectors. Each read includes the existing left edit-distance and
right candidate/extension overlap, while emission is restricted to the core
range. Results are stored at stable work indices and committed in reference
order: all plus-strand chunks for a chromosome, then all minus-strand chunks.
The former whole-chromosome scheduler remains an internal parity and benchmark
reference.

### 2. Load symbolic prefix paths

`load_prefix_hash_scan_paths` loads exact precomputed Cas9 or Cas12a paths for
the requested distance. Cas9 paths are also available for experimental prefix
lengths up to 16. If an asset is unavailable, the same paths are generated with
`build_PathTemplates`, restricted to the requested prefix and distance, then
deduplicated.

At p16, each motif has 1, 129, 7,873, 302,337, and 8,196,801 distinct
symbolic paths for d0 through d4 respectively. Paths are shared by every guide
in the query.

### 3. Build concrete hashes for each guide

Cas9 extends in the 5-prime direction, so guides are first oriented to match the
prefixHashDB template convention. Each guide is converted to a compact two-bit
alphabet.

For every symbolic path, `fill_prefix_hashes_columnwise!` selects 16 positions
from the formatted guide and folds them into a `UInt32`:

```text
hash = 0
for symbolic_position in path
    hash = (hash << 2) | guide_base[symbolic_position]
end
```

Hashes are sorted and deduplicated per guide. The 61-guide human experiment
produced about 7.0 million guide/hash associations after per-guide
deduplication.

Per-guide formatting, folding, sorting, and deduplication run as bounded tasks,
using at most `scan_threads` workers. Each task writes a distinct guide-list slot.
The heap merge remains serial and deterministic. `query_build_backend=:serial`
keeps the reference implementation; `:auto` uses it for one guide or one worker
and otherwise selects parallel construction.

### 4. Merge hashes into the compact query directory

For at most 64 guides, one bit in a `UInt64` identifies each guide. Equal hashes
from different guides are merged:

```text
concrete 16-mer hash -> 64-bit mask of compatible guides
```

The optimized lookup is not a Julia `Dict`. Sorted 32-bit hashes are split into:

- a direct bucket-offset array selected by the high hash bits;
- compact `UInt16` suffixes inside each bucket;
- a parallel `UInt64` guide-mask array.

A 26-bit presence bitmap is checked first. Most genome hashes fail this cheap
test and never access the larger directory. The bitmap may create extra work,
but cannot reject a hash that exists in the directory.

### 5. Stream FASTA chunks

The `.fai` file supplies chromosome lengths, byte offsets, and FASTA line
geometry. Each worker opens its own read-only FASTA handle and reuses one byte
buffer, one plus/minus `PrefixHashScanHit` pair, and compact candidate/radix scratch buffers.
The default logical chunk size is 2 MiB.

Each chunk includes a small overlap for the 23-base Cas9 window and the
requested edit-distance extension. FASTA newlines are removed in place in the
worker buffer. The entire chromosome is not converted to `LongDNA`.

Global `(chromosome, chunk)` work items are claimed through an atomic counter.
Stable result slots preserve chromosome, strand, and coordinate order even though
chunks finish out of order.

### 6. Detect Cas9 windows with SIMD

`scan_cas9_prefix_hits_raw_range!` first clears the worker hit vectors, then
profiles raw ASCII reference bytes in blocks. Two 32-byte AVX2 loads produce
64-bit masks describing where `A`, `C`, `G`, and `T` occur. Adjacent profiles
form a 128-base view, sufficient to evaluate 64 candidate starts together. The
allocating `scan_cas9_prefix_hits_raw_range` wrapper remains the parity reference.

Bit operations then calculate:

- starts whose complete 23-base window is unambiguous DNA;
- forward-strand windows ending in `GG`;
- reverse-strand windows beginning in `CC`, the reverse complement of `GG`.

Only set PAM bits are visited. BMI2 `PDEP` packs the corresponding 16 bases into
a `UInt32`. Reverse-strand hashes are reversed and complemented with bit
operations. This avoids allocating sequences and avoids calculating an
alignment at every genome position.

### 7. Apply the symbolic-prefix filter

Presence-bitmap checks remain in genome order. Survivors are packed as a
`UInt64` containing the 32-bit hash and local candidate start. Three stable,
worker-scratch radix passes order the full hash by its 10-bit suffix and two
11-bit bucket digits. Directory lookup then walks hashes in bucket order and
reuses the mask for repeated hashes.

Accepted hits are sorted back by candidate start before verification, preserving
the reference chromosome/strand/coordinate order. A hit produces:

```text
(candidate start, mask of potentially matching guides)
```

On the 61-guide human data, 1,547,796 genome windows passed the exact directory
filter and expanded to 1,583,279 guide/window pairs. The ratio is only 1.023
guide pairs per hit, so most hits concern one guide. `lookup_variant=:inline`
retains the former genome-order directory lookup for parity and benchmarking.
Both paths finish scanning each chunk before Myers verification and reuse all
worker buffer capacity.

### 8. Verify full edit distance with Myers

Every candidate bit is verified with
`prefix_hash_scan_raw_myers_distance`. Each guide has a precomputed Myers
equality profile. The verifier reads directly from the reusable raw FASTA
buffer, handles both strands, evaluates the 20-base guide against the extended
reference, and returns a value above the requested distance for rejected
candidates.

This verifier is allocation-free and does not construct an alignment. It is a
full Levenshtein/edit-distance filter, including indels.

### 9. Trace back accepted candidates

Only candidates whose Myers distance is at most the requested threshold are
materialized as a
`LongDNA` sequence. `align` then computes the exact alignment strings and
distance used in the output. This preserves prefixHashDB-compatible reporting
without paying traceback cost for rejected hash hits.

### 10. Commit results deterministically

Workers retain accepted hits by chromosome and strand. The main task commits
them in chromosome, strand, and position order, deduplicates complete output
records, applies early-stopping counters, and writes CSV rows.

Early stopping currently affects committed output, but workers have already
scanned and verified the genome. It therefore does not save most computation in
the streaming backend.

## Simplified pseudocode

```text
geometry = resolve_geometry(motif, distance=k, prefix=16)
paths = load_precomputed_symbolic_paths(geometry, distance=k, prefix=16)

for guide in guides:
    hashes[guide] = unique(sort(apply_each_path(paths, orient(guide))))

query = compact_directory(merge_hashes_into_guide_masks(hashes))
query = add_presence_bitmap(query, bits=26)
myers_profiles = build_myers_profiles(guides)

parallel workers claim globally scheduled overlapped FASTA chunks:
    plus_hits, minus_hits = reusable_vectors()
    plus_candidates, minus_candidates, radix_scratch = reusable_vectors()

    bases = read_and_remove_newlines(chunk)
    pam_masks = if geometry == Cas9:
        simd_find_unambiguous_NGG_and_CCN_windows(bases)
    else:
        simd_find_unambiguous_TTTV_and_BAAA_windows(bases)
    clear(hits, candidates)

    for candidate_start in set_bits(pam_masks):
        hash = bmi2_pack_oriented_16mer(bases, candidate_start)
        if presence_bitmap_contains(query, hash):
            append_packed_candidate_by_strand(
                plus_candidates, minus_candidates, hash, candidate_start)

    for strand_candidates in plus_candidates then minus_candidates:
        for hash, candidate_start in radix_order(strand_candidates):
            guide_mask = directory_lookup_or_reuse(query, hash)
            if guide_mask != 0:
                append_hit(strand_hits, candidate_start, guide_mask)
        sort_hits_by_candidate_start(strand_hits)

        for hit in strand_hits:
            for guide in set_bits(hit.guide_mask):
                if raw_myers_distance(guide, bases, hit.candidate_start) <= k:
                    alignment = traceback(
                        guide, materialize_candidate(bases, hit.candidate_start))
                    retain(alignment)

commit_retained_hits_in_reference_order()
```

## Comparison with the general path

| Stage | Optimized Cas9/Cas12a d0-d4 paths | General `:legacy` path |
|---|---|---|
| Supported query | Canonical Cas9 or Cas12a, d0-d4, 16-base hash, <=64 guides | Other motifs, hash lengths, and guide counts |
| Query structure | Presence bitmap + compact sorted directory + `UInt64` guide masks | `Dict` from hash to guide mask or guide-index vectors |
| Reference access | FAI range reads into reusable raw buffers | FASTA/2bit records loaded and converted to `LongDNA` |
| PAM search | Dedicated Cas9 or Cas12a AVX2 masks evaluate 64 starts | `findguides` over materialized chromosome sequences |
| Prefix extraction | Geometry-specific BMI2 packing from raw bytes | Sequence slicing/orientation or direct scalar hashing |
| Temporary objects | Reused raw/hit buffers and compact verified hits | More sequence objects and generic candidate ranges |
| Distance rejection | Allocation-free raw Myers before materialization | Usually `align`; some fused modes can use distance-first verification |
| Traceback | Accepted candidates only | Historically performed for many more candidate pairs |
| Parallelism | Dynamic chromosome workers | Record iteration plus backend-specific range tasks |
| Main advantage | Specialized sequential scan with cheap SIMD filtering | Generality and compatibility |

There is also an intermediate `:fused_directory` path. It uses the
compact query directory and geometry-specific direct prefix hashing, but
operates on converted chromosome sequences rather than streamed raw FASTA SIMD
blocks. It is useful as a portable fallback and correctness reference for the
fastest backend.

## Why it can beat prefixHashDB

`prefixHashDB` spends no time scanning the full genome at search time, but it
must load and traverse a large genome-derived index with many partitions and
locations. `prefixHashScan` instead performs a predictable sequential reference
pass and keeps the smaller guide-derived query structure in memory. On this
machine and workload, sequential scanning plus SIMD filtering is cheaper than
the prefixHashDB search-time index access.

This does not mean an index is intrinsically slower. Results depend on storage,
cache state, guide count, index layout, and whether index construction and
storage are included. A repeated-query workload may still favor a well-designed
genome index.

## Why it differs from Sassy

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
guides, PAM density, and candidate rate. A current fair full-human benchmark
against Rust Sassy v1 and v2 has not yet been completed.

## Current measured result

### Full distance 0-4 human sweep

The July 22, 2026 sweep used GRCh38, 24 Julia threads, 61 guides per motif,
one warmup, and five timed repetitions per algorithm and distance. Algorithm
order alternated between repetitions. Search used unlimited early-stopping
thresholds. One distance-4 prefixHashDB was built per motif with one thread and
reused for all requested distances. The timed `prefixHashScan` runs used the
normal no-statistics path; a separate untimed pass collected phase counters.

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

The d4 prefixHashDB builds took 1,419.7 s for Cas9 and 634.6 s for Cas12a.
Their resulting indexes occupied 4.60 GB and 2.95 GB. These costs are excluded
from the search table but are relevant for one-shot workflows.

Raw prefixHashDB detail output is the gold standard for this benchmark. The
initial sweep recorded `parity=false` for Cas12a only because an obsolete
Cas9-only filtering layer rejected every `extends5=false` row before comparison.
Sorting the complete raw `prefixHashDB` and `prefixHashScan` CSV rows produced
identical SHA-256 hashes at every Cas12a distance. Current parity compares exact
detail-row multisets directly, including duplicate multiplicity; prefixHashDB
rows are not independently filtered or reinterpreted. The July 23 parity
repair reports PASS with zero scan-only and zero prefix-only rows for all ten motif/distance cases.

The phase counters explain the d4 step change. Query construction took 7.61 s
for Cas9 and 7.26 s for Cas12a. Cas9 d4 expanded 8,196,801 symbolic paths into
111,720,240 guide/hash associations; Cas12a produced 109,038,602. The final
Cas9 query is about 1.09 GB. D4 is consequently both a query-construction and
scan/verification problem, rather than a small extension of d3.

Cas12a remains output-heavy. At d3 it performed 364,581 tracebacks, and at d4
it performed 1,073,287. Count-only output can bypass alignment materialization,
string construction, detail-row deduplication, and detail CSV writes for these
accepted candidates. This is a more credible large win than further Cas9
micro-optimization.

### Matched Cas9 and Cas12a human benchmark

The July 16, 2026 comparison used GRCh38, distance 3, 8 Julia threads, 61
guides per motif, unlimited early-stopping thresholds, and existing or newly
built prefixHashDB indexes. Index construction was excluded. The Cas12a set was
sampled from 61 distributed canonical `TTTV` sites in GRCh38 so every query had
a real on-target. Each algorithm was warmed once in the same process and then
measured three times. Both scans resolved to precomputed paths, `bitmask64`
queries, and `streaming_fasta_simd`.

| Motif | Results | `prefixHashScan` median | `prefixHashDB` median | Scan speedup | Scan runs | DB runs |
|---|---:|---:|---:|---:|---|---|
| Cas9 | 25,826 | 2.694 s | 22.861 s | 8.49x | 2.826, 2.694, 2.673 s | 23.263, 22.818, 22.861 s |
| Cas12a | 364,581 | 5.083 s | 14.856 s | 2.92x | 5.377, 5.083, 4.863 s | 15.114, 14.856, 14.583 s |

A separate single-pass harness measured Cas9 at 4.507 s versus 46.863 s
(10.40x) and Cas12a at 7.921 s versus 19.823 s (2.50x). These first-pass
numbers include more loading and runtime noise; the warmed medians above are
the primary comparison. All comparisons exclude prefixHashDB construction.

Both motifs had exact core-result parity on guide, distance, chromosome,
position, and strand: zero scan-only and zero prefixHashDB-only rows. The
Cas12a implementation also passed its focused scalar/SIMD/backend tests, CLI
tests, sample prefixHashDB parity benchmark, and the complete Julia test suite.

The workloads differ substantially despite equal guide counts:

| Counter | Cas9 | Cas12a |
|---|---:|---:|
| GRCh38 motif candidates | 304,418,266 | 136,277,556 |
| Exact directory hits | 1,547,796 | 1,546,515 |
| Guide/window verification pairs | 1,583,279 | 1,887,411 |
| Tracebacks and emitted rows | 25,826 | 364,581 |
| Precomputed path rows | 302,337 | 302,337 |
| Concrete query hashes before cross-guide merge | 7,044,938 | 6,947,869 |

Cas12a performs only 19% more guide/window verifications but 14.1x more
tracebacks and output commits on this guide set. Its lower speedup relative to
Cas9 is therefore primarily a verification-success and result-materialization
effect, not failure of the specialized PAM scan. Sampled CPU profiles support
this: Cas9 is dominated by the streaming scan/lookup kernel, while Cas12a has a
large additional contribution from `evaluate_prefix_hash_scan_hits!`, `align`,
sequence/string materialization, deduplication, and CSV commit. Thread/task
utilization was about 69% for Cas9 and 71% for Cas12a, leaving scheduling or
tail-latency headroom in both.

The profile's `align_ns`, `verify_ns`, and other worker fields are summed across
threads and can exceed wall time; they establish attribution, not serial stage
duration. Sampling instrumentation also raised observed wall time to 5.816 s
for Cas9 and 8.853 s for Cas12a, so those profiled timings are not used in the
speedup table.

The reusable human prefixHashDB indexes occupy approximately 3.3 GiB for Cas9
and 1.5 GiB for Cas12a. `prefixHashScan` requires neither index; it needs only
the reference FASTA and its small `.fai`.

### Earlier Cas9 scaling and tuning record

Human GRCh38, 61 Cas9 guides, distance 3, warm-cache search, pinned physical
CPUs:

| Configuration | 12 CPUs | 24 CPUs | Notes |
|---|---:|---:|---|
| Current: bucketed lookup, parallel query | 1.710 s | 1.067 s | Paired lookup A/B; no statistics |
| Inline-lookup reference, parallel query | 1.833 s | 1.139 s | Same A/B and exact output |
| Previous production confirmation | 1.725 s | 1.200 s | Before bucketed lookup |
| 2 MiB chunks, serial-query reference | 2.064 s | 1.591 s | Historical query A/B |
| 8 MiB chunks, serial-query reference | 2.313 s | 1.717 s | Chunk-confirmation run |
| `prefixHashDB` | 27.237 s | Not measured | Existing index; build excluded |

The true no-statistics path was 1.09% faster than statistics-enabled execution.
The initial end-to-end fusion comparison was 1.57% slower, but a later 12-run
prepared-query scan comparison measured 1.610 s fused versus 1.616 s buffered,
a 0.37% difference inside run variance. Fusion therefore has no demonstrated
latency effect. It reduced full-search allocated bytes from 1,113,346,392 to
1,049,027,928, a 5.78% reduction.

Reusing buffered hit vectors was tested over 15 alternating prepared-query scan
pairs. Median time improved from 1.686 s to 1.675 s (0.69%), with a paired median
delta of -33.7 ms. Allocated bytes fell from 812,254,224 to 750,851,360, a 7.56%
reduction, and every run produced the same 25,826 verified hits. Reuse passed the
memory-win/no-slowdown gate and is now the default buffered implementation. A
three-run end-to-end sanity check had a 2.579 s median, within the observed
machine variance.

Global 8 MiB chunk scheduling was tested over 15 alternating prepared-query
pairs against the former chromosome scheduler. At 12 pinned physical cores,
median scan time improved from 1.593 s to 1.407 s (13.21%), with a paired median
delta of -186.0 ms and 6.24% fewer allocated bytes. At 24 pinned physical cores,
median scan time improved from 1.417 s to 0.758 s (46.54% lower latency, or an
87.05% throughput speedup), with a paired median delta of -654.6 ms and 12.62%
fewer allocated bytes. Both modes produced the same 25,826 verified hits. The
experiment passed the 24-core 20%-improvement and 12-core no-regression gates,
so global chunk scheduling is now the production streaming scheduler. A
three-run 12-core end-to-end sanity check produced a 2.271 s median and
byte-identical 25,826-row output.

At that checkpoint, the scanner was 3.03x faster than its previous fused backend
and 10.98x faster than prefixHashDB search for this experiment. A previous
run observed peak process RSS of about 1.52 GiB. Before global scheduling, a 24-core
end-to-end median of 2.933 s was slower than the 12-core result. The scheduler
A/B result confirms whole-chromosome tail imbalance was a major cause.

These numbers are evidence for this exact machine and workload, not a universal
speed claim.

The July 2026 stabilization refactor preserved full-GRCh38 output bytes,
prepared-result signatures, semantic counters, and all 25,826 emitted rows at
both 12 and 24 CPUs. The new three-argument public API also produced identical
bytes for all 61 guides. Focused tests, the complete Julia suite, and the
Documenter build passed after separating the source files and introducing
`PrefixScanGeometry`. Refactor timing runs were made while unrelated R jobs
occupied host CPUs and were about 11% above the idle-host results in the table;
they are retained only as parity evidence and do not replace the uncontaminated
production measurements.

A staged GRCh38 sweep varied chunk size (2, 4, 8, and 16 MiB), prefilter
width (22, 24, and 26 bits), and bucket prefix (9 through 12 bases) at 12 and 24
pinned physical CPUs. Every configuration preserved result signatures, output
bytes, and semantic statistics counters. The initial sweep identified 2 MiB /
26 bits / 11 bases as a memory-win candidate.

A follow-up used 15 alternating full-GRCh38 pairs and five allocation runs per
configuration. At 12 CPUs, 2 MiB reduced median end-to-end time from 2.313 s to
2.227 s (3.73%) and cumulative allocation from 1.010 GB to 0.869 GB (14.0%). At
24 CPUs, time fell from 1.717 s to 1.682 s (2.00%) and allocation from 1.153 GB
to 0.895 GB (22.3%). Three isolated 24-CPU processes per configuration measured
median peak RSS of 1.650 GB for 8 MiB and 1.393 GB for 2 MiB, a 15.6% reduction.
All runs produced the same 25,826 results, bytes, and semantic counters. The
2 MiB configuration passed the memory-win/no-slowdown gate and is now the
production default.

Serial query construction was then compared with bounded per-guide tasks while
retaining the deterministic heap merge. Over 15 alternating 61-guide GRCh38
pairs, 12-CPU median query construction improved from 0.728 s to 0.320 s (56.1%)
and end-to-end time from 2.064 s to 1.725 s (16.4%). At 24 CPUs, query
construction improved from 0.778 s to 0.343 s (56.0%) and end-to-end time from
1.591 s to 1.200 s (24.6%). Allocations changed by less than 0.01%, and all
result signatures, output bytes, and semantic counters matched.

The 8-guide crossover also improved end-to-end time by 5.6% at 12 CPUs and 8.7%
at 24 CPUs. One-guide `:auto` execution stays serial. Parallel per-guide query
construction therefore passed its promotion gate and is now the `:auto` behavior
for multi-guide compact queries. Reproduce it with:

```bash
CHOPOFF_TUNING_STAGE=query \
  julia --project=. scripts/benchmark_prefix_hash_scan_tuning.jl
```

The promoted lookup experiment used 15 alternating full-GRCh38 pairs and five
allocation runs per variant. At 12 CPUs, radix/bucket ordering reduced prepared
scan time from 1.321 s to 1.171 s (11.3%) and end-to-end time from 1.833 s to
1.710 s (6.7%). Allocated bytes rose from 0.869 GB to 0.885 GB (1.8%), while
three isolated processes reduced median peak RSS from 1.092 GB to 1.063 GB
(2.6%). At 24 CPUs, prepared scan time fell from 0.717 s to 0.649 s (9.5%) and
end-to-end time from 1.139 s to 1.067 s (6.2%). Allocated bytes rose from 0.896
GB to 0.926 GB (3.4%), while median peak RSS fell from 1.362 GB to 1.259 GB
(7.6%). Every run produced the same 25,826 results, output bytes, and semantic
counters. The experiment passed its 5% latency and 5% memory-regression gates,
so compatible `:auto` searches now use bucketed lookup. Reproduce it with:

```bash
CHOPOFF_TUNING_STAGE=lookup \
  julia --project=. scripts/benchmark_prefix_hash_scan_tuning.jl
```

The current query contains 6,898,183 exact hashes. Its offset, suffix/mask, and
presence arrays occupy about 16.8 MB, 69.0 MB, and 8.4 MB respectively. The
26-bit bitmap has 3,121,043 set prefixes (4.65%), implying roughly 14.2 million
prefilter survivors from the 304.4 million GRCh38 PAM candidates before exact
directory lookup.

Pre-experiment 2 MiB/parallel-query production sampling attributed 35.8%
cumulative samples to SIMD scan/lookup and 19.5% to directory lookup at 12
CPUs. At 24 CPUs, those shares fell to 16.0% and 7.9%. Worker wait rose from
24.2% to 44.4%. Myers verification accounted for only 2.6% and 1.8%; FASTA
range reading accounted for 3.2% and 1.5%.
Hardware performance counters remain unavailable, so these are Julia sampling
shares rather than direct cache-miss measurements.

Earlier statistics-enabled profiling attributed about 0.97 s of a 2.51 s
12-core run to serial query construction and 1.54 s to scan plus commit. The new
parallel result removes query construction as the primary bottleneck; prepared
scanning now dominates. Hardware performance counters remain unavailable on
this host.

## Current limitations and issues

1. The fastest kernels are separate Cas9 and Cas12a geometries at distances 0
   through 4 with a 16-base prefix; other motifs and prefix lengths use slower
   paths.
   Shared validation, bounds, overlap, and scheduling use
   `PrefixScanGeometry{Kind}` without moving motif branches into SIMD loops.
2. A `UInt64` guide mask limits the fused path to 64 guides.
3. AVX2 and BMI2 are required; there is no equivalent ARM or AVX-512 kernel.
4. FASTA requires `.fai`; optimized streaming is not implemented for 2bit.
5. Ambiguous query guides are rejected.
6. The supported API intentionally skips any candidate with a non-ACGT base in
   its complete 23-base Cas9 or 25-base Cas12a guide/PAM window. There is no
   ambiguous-reference fallback.
7. Early stopping does not cancel already scheduled scan/verification work.
8. Global scheduling adds per-chunk result containers and concurrent FASTA
   seeks. Current evidence is for warm-cache GRCh38; cold-cache and networked
   filesystems have not been measured.
9. Buffered workers retain the largest observed per-chunk hit and prefilter-survivor
   capacities. The extra candidate/radix buffers improve locality but explain the
   1.8-3.4% increase in cumulative allocation.
10. Requested statistics still add counters and `time_ns()` calls to hot loops;
    `stats=nothing` now compiles those operations out.
11. Query construction is rebuilt for every call. Per-guide hash lists are now
    parallel, but their heap merge and directory construction remain serial;
    cross-run reuse is intentionally out of scope for the one-shot workload.
12. Compact-directory accesses are now bucket ordered, but the 8.4 MB presence
    bitmap is still probed in genome order and may remain memory-latency bound.
13. `PrefixHashScanStats` contains summed worker CPU times and wall-clock fields;
    those values cannot be compared directly without clear labeling.
14. The exported three-argument `search_prefixHashScan` and standalone CLI
    support Cas9 and Cas12a at distances 0 through 4. Advanced tuning remains
    on the experimental
    four-argument method.
15. Source separates constant-free helpers, Cas9 and Cas12a kernels,
    verification, streaming, and orchestration in the parent `CHOPOFF` module.

## Potential speed optimizations

### Highest-priority experiments

1. **Completed: profile the current production path.** At 12 CPUs, SIMD
   scan/lookup and directory lookup accounted for 35.8% and 19.5% cumulative
   samples. At 24 CPUs these fell to 16.0% and 7.9%, while wait rose to 44.4%.
   Hardware counters were unavailable.
2. **Completed: true no-statistics hot path.** Per-worker statistics, counters,
   and timers are absent when `stats=nothing`; measured gain was 1.09%.
3. **Completed experiment: fuse lookup and verification.** Immediate
   verification removed `PrefixHashScanHit` vectors and 5.78% of allocated bytes,
   but produced no measurable latency change. The buffered implementation remains
   the `:auto` backend.
4. **Completed: reuse buffered hit vectors.** Worker-owned scratch vectors cut
   prepared-scan allocations by 7.56% and improved median scan time by 0.69%
   without changing output.
5. **Completed: schedule chunks globally.** Stable `(chromosome, chunk)` work
   items improved prepared scan time by 13.21% at 12 cores and reduced 24-core
   latency by 46.54%. Exact ordering, statistics counters, and output parity are
   preserved; this is now the production scheduler.

### Cache and lookup experiments

6. **Completed: sweep `prefilter_bits`, `bucket_bases`, and chunk size.** The
   staged GRCh38 sweep retained 26 bits and 11 bases while identifying 2 MiB as
   the memory candidate. Reproduce with
   `scripts/benchmark_prefix_hash_scan_tuning.jl` and
   `CHOPOFF_TUNING_STAGE=chunk|prefilter|bucket|final`.
7. **Completed: confirm and promote 2 MiB chunks.** Fifteen paired runs at both
   CPU counts showed 2.0-3.7% lower latency, 14.0-22.3% less cumulative
   allocation, and 15.6% lower median peak RSS without changing output. Two MiB
   is now the production default.
8. **Completed: radix-order directory lookup.** Worker-scratch radix passes
   over prefilter survivors improved prepared scans by 9.5-11.3% and end-to-end
   time by 6.2-6.7%. Allocation rose 1.8-3.4%, peak RSS fell 2.6-7.6%, and exact
   parity held. Compatible `:auto` searches now use it.
9. Compare the current bitmap/directory with a cache-conscious blocked Bloom,
   xor, or quotient-style prefilter. Any replacement must be exact after the
   final directory lookup and must be evaluated at equal memory usage.
10. Measure NUMA placement and pinning. Per-socket query replicas may outperform
    shared random access if the directory is frequently fetched across sockets.

### Query and verification experiments

11. **Completed: parallelize per-guide concrete hash lists.** Bounded tasks
    retain the deterministic serial merge. Query construction improved about 56%
    and 61-guide end-to-end time improved 16.4-24.6%; `:auto` now selects it for
    multi-guide compact queries.
12. SIMD-vectorize raw Myers or traceback across independent candidate hits.
    Cas9 averages only 1.023 guide pairs per hash hit and remains scan-heavy,
    but the Cas12a human workload emitted 364,581 accepted candidates and made
    verification/materialization substantial. Benchmark the motifs separately;
    a Cas12a win must not complicate or regress the Cas9 path.
13. Add an AVX-512 profiling kernel that handles more input bytes per iteration.
    Expect an incremental gain unless profiling proves the SIMD base-profile
    stage dominates.

### Lower-priority micro-optimizations

14. Reduce FASTA newline-compaction copies or scan an mmap-backed representation
    directly. This adds format complexity and should follow I/O measurements.
15. Avoid string conversion and generic dedup keys during final commit. This is
    low priority for the 25,826-row Cas9 workload but credible for the
    364,581-row Cas12a workload.

## Usability and feature work

Performance work should not be mixed blindly with generalization. Useful next
features are:

1. **Completed:** export and document Cas9 and Cas12a distances 0 through 4 for
   `search_prefixHashScan`, add a standalone CLI search mode, and report the
   resolved backend and execution modes without enabling statistics.
2. Support more than 64 guides in one genome pass with a compact hash-to-guide
   ID list. Rescanning the genome in 64-guide batches would be simpler but loses
   the main amortization advantage.
3. Make early stopping cancel future chunks while preserving deterministic
   output semantics.
4. **Completed:** add distance 4 at p16 as the functional upper bound and
   benchmark mode.
5. **Completed:** add Cas12a as a separate specialized scan geometry using its
   precomputed paths without forcing Cas9 constants into a generic SIMD loop.
6. Add an optimized 2bit reader and define exact IUPAC-reference behavior.

## Recommended next work after the d0-d4 sweep

The following order balances correctness, user-visible value, and measured
performance. Performance changes retain the existing 10% GRCh38 improvement
gate unless they provide an independently valuable feature.

1. **Completed: use prefixHashDB as the benchmark oracle.** Compare raw exact
   detail-row multisets, including duplicate multiplicity, without an
   intermediate motif-specific validator. Keep core-locus comparison only for
   external Sassy implementations whose traceback strings may differ.
2. **Implement count-only output.** Add an `output=:detail|:counts` API and CLI
   surface in `src/db_prefix_hash_scan.jl`, aggregate raw Myers distances in
   `src/prefix_hash_scan/verification.jl`, and bypass traceback, strings,
   detail deduplication, and detail CSV emission. Tests in
   `test/src/prefix_hash_scan.jl` must compare complete counts with
   `summarize_offtargets(detail)` across motifs, distances, strands, indels, and
   chunk boundaries. Capped searches must report `complete=false`. Cas12a d3
   and d4 are the primary performance gates.
3. **Run a d4 prefix-length and memory sweep before redesigning the query.**
   Compare p14, p15, and p16 on GRCh38 and record query-build time, peak RSS,
   unique hashes, directory hits, verification calls, and end-to-end latency.
   Exact detail parity is mandatory. Promote a different d4 prefix only if it
   reduces latency or memory by at least 10%; otherwise retain functional p16.
4. **Optimize Cas12a detail commit if count mode confirms the attribution.**
   Target `evaluate_prefix_hash_scan_hits!`,
   `commit_prefix_hash_scan_verified!`, compact deduplication keys, delayed
   string construction, and possibly batched traceback. Measure Cas9 separately
   and allow no more than 3% Cas9 d3 regression.
5. **Support more than 64 guides in one reference pass.** Preserve the current
   `UInt64` fast path and add `sorted_hashes`, `hash_offsets`, and
   `flat_guide_ids` for larger sets. Demonstrate bounded construction and exact,
   deterministic output for 64, 256, 1,024, and 4,096 guides without rescanning
   the genome in 64-guide batches.
6. **Make early stopping cancel unclaimed chunks.** Preserve deterministic
   detail output and mark partial count output incomplete. This work is valuable
   only for capped searches and should not complicate the complete-search hot
   path.
7. **Keep hardware-specific work evidence-driven.** NUMA pinning and per-socket
   query replicas are the first Cas9 experiments; cache-conscious prefilters,
   AVX-512, mmap/2bit streaming, and ARM SIMD remain follow-ups. Do not replace
   the current path without exact parity and a credible end-to-end gain.

The simplest high-value sequence is therefore: retain direct prefixHashDB
parity, ship count mode, then use the d4 prefix sweep to decide whether
compressed or staged d4 query research is justified.

## Recommended decision rule

Continue speed research while an experiment has a credible route to at least a
10% end-to-end improvement on GRCh38. No-statistics execution, scan/verify
fusion, buffered-vector reuse, global scheduling, parameter sweeps, 2 MiB
chunks, parallel guide hashing, production profiling, and radix-ordered lookup
are now measured.

Cache-local directory lookup passed, but its paired end-to-end gain was 6.2-6.7%
and lookup is only 7.9% of cumulative samples at 24 CPUs. Cas12a specialization
is now complete and confirms that motif geometries should remain separate.
Future Cas9 work should target genome-order presence-bitmap probes or explain
worker wait through NUMA/pinning measurements. Future Cas12a work should first
target accepted-candidate verification, alignment materialization, deduplication,
and commit. Require exact parity and report both motifs so an optimization for
the result-heavy Cas12a workload does not regress Cas9.

The likely remaining gain without a genome index is meaningful but smaller than
the previous 10x improvement. Another 1.5-3x is plausible only if profiling
confirms avoidable lookup, scheduling, or temporary-data costs. Another 10x is
unlikely because every exact indexless search must still inspect the reference.

## Open questions

- How much of the remaining 0.32-0.34 s query build is the serial heap merge
  and directory construction?
- Are the remaining genome-order presence-bitmap probes limited by last-level
  cache misses or memory latency?
- How much of the 24-core 44.4% wait share is tail imbalance versus task/runtime
  overhead, and can NUMA-local query placement reduce it?
- Why does immediate verification show no latency gain despite allocating 5.78%
  fewer bytes: instruction pressure or phase-locality loss?
- How much Cas12a time can be removed by batching traceback or replacing
  string-based commit/deduplication while preserving exact output ordering?
- Why are Cas12a query/path preparation costs higher than Cas9 despite nearly
  identical concrete-hash and path counts?
- Does Cas12a retain its 2.9x advantage for guide sets with fewer accepted
  off-targets, where traceback and CSV output do not dominate?

## Generalization roadmap: primary indexless search

This July 18, 2026 roadmap makes `prefixHashScan` the intended primary CHOPOFF
algorithm for ordinary reference-genome search. The optimization record above
is retained as historical evidence; the priorities below supersede its
speed-first ordering.

The target is an exact production search with:

- no genome-specific CHOPOFF database or build step;
- only reusable motif path assets and the standard FASTA `.fai`;
- one reference scan for the complete guide set;
- deterministic results with prefixHashDB-compatible coordinates and distances;
- specialized hot loops for important motif geometries, without motif branches
  inside those loops.

The next work should generalize the useful surface before attempting more
Cas9 micro-optimization. Count output is the next implementation target.

Implementation order changed July 22, 2026: detailed output for distances 0
through 3 shipped first at p16, followed by distance 4 as the functional and
benchmarking ceiling.

### Priority 1: production distances 0 through 4 (completed at p16)

Distances 0, 1, 2, 3, and 4 are supported public configurations for Cas9 and
Cas12a. A lower requested distance does not run a larger-threshold query and
discard higher-distance results. Distances 0 through 3 load exact single-file
p16 assets; distance 4 loads the existing canonical split matrices.

Distance parameterizes path selection, chunk overlap, and the Myers threshold.
The existing Cas9 and Cas12a PAM/profile kernels remain separate and are not
copied once per distance. Dispatch resolves the distance before entering each
hot loop. Production prefix length remains fixed at 16.

Future prefix-length tuning should use a GRCh38 sweep of query-build time, peak
memory, exact directory hits, verification calls, and end-to-end latency.
Exact detail parity with prefixHashDB is required at every supported distance;
count/detail parity becomes a release gate when count mode is implemented.

### Priority 2: count output mode (next)

Only two output modes are needed:

1. `detail` retains the current output with guide, alignment strings, exact
   distance, chromosome, coordinate, and strand. Accepted candidates require
   traceback.
2. `counts` returns one summary row per guide with `D0` through `Dk` and a
   `complete` flag. It does not produce location-only rows.

Count mode should use the raw Myers result after the prefix filter and avoid
candidate materialization, alignment strings, traceback, detail-row
deduplication keys, and detail CSV writes whenever correctness permits. Its
counts must match `summarize_offtargets(detail; distance=k)` for a complete
search. This parity is a release gate rather than an assumed property: it must
cover indels, both strands, chunk boundaries, and sites with competing optimal
alignments.

`complete=true` means the whole reference was processed for that guide.
Explicit early-stopping limits may return capped counts, but then `complete`
must be false. A count must never be presented as exact after work for that
guide was cancelled. Early stopping should eventually cancel unclaimed chunks,
instead of only suppressing records during the deterministic commit.

This mode is also a performance feature. In the current 61-guide human
experiments, detail output performs 25,826 Cas9 tracebacks and 364,581 Cas12a
tracebacks. Count mode should bypass both workloads while preserving exact
per-distance counts.

### Priority 3: thousands of guides in one pass

The current `UInt64` guide mask should remain the fast representation for at
most 64 guides. It should not be widened into a dense mask whose size is paid
for at every hash, and the genome should not be rescanned in 64-guide batches.

For larger queries, use the same sorted hash directory with a compact
hash-to-guide-ID representation:

```text
sorted_hashes
hash_offsets
flat_guide_ids
```

The guide-ID element width can be selected from the guide count. Singleton
hashes may store one guide ID inline if measurement justifies the added
representation. Directory lookup should return a contiguous guide-ID slice,
which the verifier consumes without allocating.

Construction must remain bounded and parallel. Per-guide hash lists can be
built independently, but merging should stream into the compact directory
rather than retain unnecessary duplicate associations. The main scale series
is 64, 256, 1,024, and 4,096 guides, recording:

- query-build wall time and peak RSS;
- unique hashes and guide/hash associations;
- scan, lookup, verification, and output time;
- exact-result parity and deterministic ordering;
- warm-cache and cold-cache end-to-end latency.

The architecture succeeds only if each guide set is served by one genome pass.
If the concrete path expansion becomes the limiting memory cost, reduce or
stage the query representation rather than silently batching reference scans.

### Priority 4: bounded IUPAC ambiguity

Ambiguity support is required, but simple expansion is only one part of the
solution. Query, prefix, PAM, and verification ambiguity have different costs:

- An ambiguous query base can expand to its compatible A/C/G/T choices while
  constructing concrete prefix hashes.
- An ambiguous reference base inside the hashed prefix can generate compatible
  concrete hashes before directory lookup.
- PAM ambiguity should be evaluated directly with IUPAC compatibility masks,
  not by materializing every PAM string.
- Ambiguity outside the hashed prefix requires an IUPAC-aware Myers profile and
  exact final verification.

All expansion must be bounded by `motif.ambig_max` and use CHOPOFF's existing
`iscompatible` semantics. A candidate containing `a` independently ambiguous
bases can require up to `4^a` concrete hashes, so the optimized path needs an
explicit expansion bound and a correct scalar/fallback path above it. It must
never reject a valid candidate merely because expansion would be expensive.

Implement ambiguity first in the portable correctness backend, establish
linearDB/prefixHashDB parity, and specialize SIMD handling only after candidate
rates and ambiguity distributions are measured.

### Completed: distance 4 benchmark ceiling

Distance 4 is supported for Cas9 and Cas12a at p16 as the largest functional
search threshold and as a benchmarking opportunity. It uses the existing split
d4 matrices: 8,196,801 symbolic paths and about 125 MiB per motif. The loader
concatenates those canonical matrices directly and skips the generic deduplication
and threshold-filtering work.

This mode intentionally reuses the d0-d4 compact query and motif-specific scan
architecture. It is not expected to match lower-distance speed or memory use. In
an earlier isolated 61-guide Cas9 memory run with 24 query workers, query
construction took 10.8 seconds, produced 111,720,240 guide/hash associations
and a 1.09 GB final query, reaching 4.51 GB peak process RSS. No guide-count cap below the existing 64 is imposed; callers must provision memory accordingly. Future compressed or
staged d4 representations remain benchmark research, not a release requirement.

### Workload-specific performance after generalization

Cas9 and Cas12a now have different remaining costs and should have separate
performance gates:

- Cas9 work should investigate genome-order presence-bitmap latency,
  last-level-cache behavior, NUMA-local query copies, and scheduler tail
  imbalance.
- Cas12a detail work should target accepted-hit materialization, traceback,
  string construction, deduplication, and CSV commit.
- Count mode should be benchmarked independently because it removes most of the
  Cas12a-specific output cost and exposes the scan/verification ceiling.
- AVX-512, ARM SIMD, 2bit streaming, and further prefilter experiments remain
  follow-up portability or speed work, not prerequisites for the first
  generalized release.

Continue an optimization only when it has a credible route to at least 10%
end-to-end improvement in its intended workload. Always report detail and
count modes separately and verify that a motif-specific win does not regress
the other production geometry.

### Replacement gates for prefixHashDB

`prefixHashScan` can replace prefixHashDB as the documented default reference
search when all of the following hold:

1. Cas9 and Cas12a distances 0-4 have exact detail parity on sample, randomized,
   and human-scale fixtures.
2. Complete count output exactly matches summaries derived from detail output;
   capped output is marked `complete=false`.
3. The one-pass large-guide representation is demonstrated through at least
   4,096 guides with bounded peak memory.
4. Scalar and optimized backends have identical results, and unsupported CPU or
   input configurations select a correct fallback rather than fail silently.
5. Warm-cache, cold-cache, thread-scaling, query-build, scan, verification, and
   output costs are reported separately.
6. Existing Cas9/d3 and Cas12a/d3 detail latency regresses by no more than 3%
   unless the change provides a measured feature-level benefit that justifies
   it.

After these gates pass, update the README, public documentation, and CLI
examples to lead with the no-build workflow. Keep prefixHashDB available during
one compatibility period for validation, repeated-query workloads, VCF-related
workflows, and configurations not yet supported by the scan path.
