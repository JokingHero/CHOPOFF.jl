# prefixHashScan

## Status and scope

`prefixHashScan` is an indexless CRISPR off-target search prototype. It reuses
the symbolic prefix paths from `prefixHashDB`, but builds a guide-specific query
structure in memory and scans the reference genome directly.

This document focuses on the current optimized configuration:

- Cas9 with an `NGG` PAM
- 20-base guides
- edit distance 3
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

## Core idea

The algorithm moves part of the alignment work to the query side:

1. Enumerate every symbolic way in which a 16-base reference prefix can be
   produced from a guide while spending at most three edits.
2. Apply those symbolic paths to each concrete guide and encode the resulting
   16-mers as 32-bit integers.
3. Scan only Cas9-compatible genome windows and test whether their 16-mer is in
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

## Optimized Cas9/d3 path

### 1. Select the backend

`search_prefixHashScan(...; scan_backend=:auto)` selects
`:streaming_fasta_simd` when all of the following hold:

- the motif has the Cas9 guide/PAM geometry;
- `distance == 3`;
- `hash_len == 16`;
- the query uses the 64-bit guide mask representation;
- there are no more than 64 guides;
- the reference is FASTA;
- AVX2 and BMI2 are available.

If Cas9/d3 is supported but raw FASTA SIMD is unavailable, `:auto` selects the
intermediate `:fused_directory` backend. Otherwise it selects `:legacy`.

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

`load_prefix_hash_scan_paths` loads precomputed Cas9 paths when available. For
Cas9, paths are available for prefix lengths up to 16 and distances up to 4.
If an asset is unavailable, the same paths are generated with
`build_PathTemplates`, restricted to the requested prefix and distance, then
deduplicated.

For Cas9/d3/hash-length 16, the current asset contains 302,337 distinct symbolic
paths. These paths are shared by every guide in the query.

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
three-base edit-distance extension. FASTA newlines are removed in place in the
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
reference, and returns a value above 3 for rejected candidates.

This verifier is allocation-free and does not construct an alignment. It is a
full Levenshtein/edit-distance filter, including indels.

### 9. Trace back accepted candidates

Only candidates whose Myers distance is at most 3 are materialized as a
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
paths = load_precomputed_symbolic_paths(Cas9, distance=3, prefix=16)

for guide in guides:
    hashes[guide] = unique(sort(apply_each_path(paths, orient(guide))))

query = compact_directory(merge_hashes_into_guide_masks(hashes))
query = add_presence_bitmap(query, bits=26)
myers_profiles = build_myers_profiles(guides)

parallel workers claim globally scheduled overlapped FASTA chunks:
    plus_hits, minus_hits = reusable_vectors()
    plus_candidates, minus_candidates, radix_scratch = reusable_vectors()

    bases = read_and_remove_newlines(chunk)
    pam_masks = simd_find_unambiguous_NGG_and_CCN_windows(bases)
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
                if raw_myers_distance(guide, bases, hit.candidate_start) <= 3:
                    alignment = traceback(
                        guide, materialize_candidate(bases, hit.candidate_start))
                    retain(alignment)

commit_retained_hits_in_reference_order()
```

## Comparison with the general path

| Stage | Optimized Cas9/d3 path | General `:legacy` path |
|---|---|---|
| Supported query | Cas9, d3, 16-base hash, <=64 guides | Other distances, motifs, hash lengths, and guide counts |
| Query structure | Presence bitmap + compact sorted directory + `UInt64` guide masks | `Dict` from hash to guide mask or guide-index vectors |
| Reference access | FAI range reads into reusable raw buffers | FASTA/2bit records loaded and converted to `LongDNA` |
| PAM search | 64 candidate starts evaluated with AVX2 bit masks | `findguides` over materialized chromosome sequences |
| Prefix extraction | BMI2 packing from raw bytes | Sequence slicing/orientation or direct scalar Cas9 hashing |
| Temporary objects | Reused raw/hit buffers and compact verified hits | More sequence objects and generic candidate ranges |
| Distance rejection | Allocation-free raw Myers before materialization | Usually `align`; some fused modes can use distance-first verification |
| Traceback | Accepted candidates only | Historically performed for many more candidate pairs |
| Parallelism | Dynamic chromosome workers | Record iteration plus backend-specific range tasks |
| Main advantage | Specialized sequential scan with cheap SIMD filtering | Generality and compatibility |

There is also an intermediate Cas9/d3 `:fused_directory` path. It uses the
compact query directory and fused rolling Cas9 scan, but operates on converted
chromosome sequences rather than streamed raw FASTA SIMD blocks. It is useful
as a portable fallback and correctness reference for the fastest backend.

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

1. The fastest backend is hard-coded around Cas9/d3/16-base-prefix geometry.
2. A `UInt64` guide mask limits the fused path to 64 guides.
3. AVX2 and BMI2 are required; there is no equivalent ARM or AVX-512 kernel.
4. FASTA requires `.fai`; optimized streaming is not implemented for 2bit.
5. Ambiguous query guides are rejected.
6. Raw SIMD scanning requires an unambiguous 23-base candidate window. This is
   tested against the existing fused behavior, but broader IUPAC policy should
   be specified explicitly before treating the backend as a public interface.
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
14. `search_prefixHashScan` is not exported, documented as public API, or
    available through the standalone CLI.
15. The implementation remains concentrated in one large source file, which
    makes kernel-level profiling and maintenance harder.

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
12. SIMD-vectorize raw Myers across independent candidate hits. The current
    average of 1.023 guides per hash hit makes SIMD across guides unattractive;
    batching hits for the same guide is more plausible. Do this only if current
    profiling shows Myers is substantial.
13. Add an AVX-512 profiling kernel that handles more input bytes per iteration.
    Expect an incremental gain unless profiling proves the SIMD base-profile
    stage dominates.

### Lower-priority micro-optimizations

14. Reduce FASTA newline-compaction copies or scan an mmap-backed representation
    directly. This adds format complexity and should follow I/O measurements.
15. Avoid string conversion and generic dedup keys during final commit. Only
    25,826 rows were emitted in the human run, so this is unlikely to matter.

## Usability and feature work

Performance work should not be mixed blindly with generalization. Useful next
features are:

1. Export and document `search_prefixHashScan`, add a standalone CLI search
   mode, and report which backend `:auto` selected.
2. Support more than 64 guides in one genome pass with a compact hash-to-guide
   ID list. Rescanning the genome in 64-guide batches would be simpler but loses
   the main amortization advantage.
3. Make early stopping cancel future chunks while preserving deterministic
   output semantics.
4. Add specialized distance 1, 2, and 4 kernels and benchmark prefix lengths
   independently for each distance.
5. Add Cas12a as a separate specialized scan geometry using its precomputed
   paths. Do not force Cas9 constants into a nominally generic SIMD kernel.
6. Add an optimized 2bit reader and define exact IUPAC-reference behavior.

## Recommended decision rule

Continue speed research while an experiment has a credible route to at least a
10% end-to-end improvement on GRCh38. No-statistics execution, scan/verify
fusion, buffered-vector reuse, global scheduling, parameter sweeps, 2 MiB
chunks, parallel guide hashing, production profiling, and radix-ordered lookup
are now measured.

Cache-local directory lookup passed, but its paired end-to-end gain was 6.2-6.7%
and lookup is only 7.9% of cumulative samples at 24 CPUs. The next speed work
should target the remaining genome-order presence-bitmap probes or explain the
44.4% worker-wait share through NUMA/pinning measurements. Require equal-memory
comparisons and exact final directory verification. If neither has a credible
10% route, stabilize Cas9/d3 before expanding motifs.

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
- What exact ambiguous-reference semantics should the supported API guarantee?
