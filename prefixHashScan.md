# prefixHashScan

## Status and scope

`prefixHashScan` is an indexless CRISPR off-target search. Its public Julia API
accepts registered motifs and custom `Motif` objects at edit distances 0 through
4. Hand-written Cas9 and Cas12a kernels remain the canonical fast paths. Other
eligible motifs use a motif-specialized generic kernel; configurations outside
that kernel's envelope use the exact legacy engine. The search reuses the
symbolic prefix paths from `prefixHashDB`, builds a guide-specific query
structure in memory, and scans the reference genome directly.

The current optimized envelope is:

- Cas9: 20-base guide with an `NGG` PAM;
- Cas12a: 21-base guide with a `TTTV` PAM;
- generic: 16 through 64 guide bases, complete motif span at most 65 bases;
- one contiguous PAM block at any position, or no PAM;
- forward, reverse, or both strands;
- edit distance 0 through 4
- 16-base prefix
- `guide length - distance >= 16`;
- 64 guides per optimized query batch; larger lists are batched automatically
- FASTA with a standard `.fai`, or a `.2bit` reference
- x86 CPU with AVX2 and BMI2; AVX-512F/BW is used on qualified CPUs

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

The public API and CLI accept larger guide lists and search them as sequential
64-guide batches. Detail mode appends batches to a sibling staging file and
atomically publishes it after every batch succeeds. Count mode merges batch
matrices and writes one row per unique guide after applying deterministic
per-distance caps. This is intentional: a GRCh38 benchmark found batching
about three times faster than the tested one-pass large-guide representations.

"Indexless" means that no CHOPOFF genome database must be built. The optimized
FASTA reader still requires the small, standard `.fai` random-access index.

The implementation is split by stable responsibility:

- `src/db_prefix_hash_scan.jl`: shared types, orchestration, and public API;
- `src/prefix_hash_scan/query.jl`: symbolic paths, hashes, directory, prefilter;
- `src/prefix_hash_scan/kernel_common.jl`: geometry-neutral SIMD/lookup primitives;
- `src/prefix_hash_scan/cas9.jl`: scalar and typed x86 SIMD Cas9 scan kernels;
- `src/prefix_hash_scan/cas12a.jl`: scalar and typed x86 SIMD Cas12a scan kernels;
- `src/prefix_hash_scan/generic.jl`: compiled motif-specialized generic kernel;
- `src/prefix_hash_scan/verification.jl`: Myers, traceback, and result commit;
- `src/prefix_hash_scan/streaming.jl`: FASTA/2bit streaming and global scheduler.

`PrefixScanGeometry{Kind,Matcher}` supplies guide, PAM, prefix, distance,
candidate-span, overlap, and optional compiled motif matching to validation and
orchestration. Literal Cas9 and Cas12a constants stay inside separate SIMD hot
loops.

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
the requested distance-specific geometry before dispatch. Scalar or x86 SIMD
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
same FASTA, guide-count, and x86 SIMD requirements as Cas9. The public API
selects it with `motif="Cas12a"`; the CLI uses `--motif Cas12a`.

### Typed generic specialization

`resolve_generic_prefix_scan_geometry` converts eligible motif properties into
the `PrefixScanMatcher` type: enabled strands, constrained IUPAC positions,
guide offsets after PAM removal, orientation, and coordinate offsets. Generated
functions turn that type into straight-line validity, PAM-matching, and prefix
packing operations. There are no runtime motif branches inside the 64-start
SIMD block loop.

Guide and motif lengths are not fixed to Cas9 or Cas12a. For example, a
25-base guide followed by an `NNT` PAM resolves to
`PrefixScanGeometry{:generic}` with a 28-base candidate span and a 16-base
prefix at distances 0 through 4:

```julia
motif = Motif(
    "25N_NNT",
    repeat("N", 25) * "NNT",
    repeat("X", 25) * "NNT",
    true, true, 4, true, 0,
)
```

An all-`X` PAM description produces an empty PAM range and therefore a PAMless
search. The optimized generic envelope requires a 16-through-64-base guide, a
complete span no longer than 65 bases, distance 0 through 4, and
`guide length - distance >= 16`. Other valid motifs remain correct through
`:legacy`.

## Optimized streaming path

### 1. Select the backend

`search_prefixHashScan(...; scan_backend=:auto)` selects
`:streaming_fasta_simd` when all of the following hold:

- the motif has a canonical or eligible typed generic geometry;
- the distance is between 0 and 4;
- `hash_len == 16`;
- the query uses the 64-bit guide mask representation;
- there are no more than 64 guides;
- the reference is FASTA;
- AVX2/BMI2 or AVX-512F/BW/BMI2 are available.

`simd_backend=:auto` selects AVX-512 only for benchmark-qualified CPU-family
and specialized-geometry pairs; generic geometries retain AVX2. `:avx2` and
`:avx512` explicitly force a supported ISA.
If raw FASTA SIMD is unavailable, `scan_backend=:auto` selects the intermediate
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

`load_prefix_hash_scan_paths` loads exact precomputed paths by guide length,
distance, and prefix length. Existing 20-base and 21-base assets are reused
regardless of PAM sequence or position. If an asset is unavailable, as for a
25-base guide, paths are generated once before guide batching with
`build_PathTemplates`, restricted to the requested prefix and distance,
deduplicated, and reused by every batch. Path generation does not select the
scan kernel and does not force the generic engine to `:legacy`.

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

#### Large-guide batching and benchmark decision

Guide lists larger than 64 are partitioned in input order. Every batch runs the
normal `UInt64` query path, and its rows are appended without another CSV
header. Early stopping remains independent per guide. Output order is
deterministic and batch-major.

The selection was measured on GRCh38 with 1,024 Cas9/d3 guides, eight threads,
seven rotated timed repetitions, and exact detail-row multiset comparison:

| workload | sequential 64-guide batches | best large directory | ratio |
|---|---:|---:|---:|
| dispersed guides | 608 s [556, 628] | 1,810 s [1,722, 1,977] | 2.98× |
| related guides | 407 s [369, 420] | 1,130 s [933, 1,299] | 2.78× |

Brackets are bootstrap 95% intervals for the median. Every comparison passed
exact output parity. Wider presence filters, a 12-base directory bucket, and a
reference-aware two-pass query reduced query memory in some cases but did not
close the end-to-end gap. The large-directory prototype was therefore removed.
It should only be reconsidered after a bounded-memory design beats batching by
at least 10% on both workloads.

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
profiles raw ASCII reference bytes in blocks. AVX2 uses two 32-byte loads;
AVX-512BW uses one 64-byte load and mask comparisons. Both produce identical
64-bit `A`, `C`, `G`, and `T` profiles. Adjacent profiles form a 128-base view,
sufficient to evaluate 64 candidate starts together. The allocating
`scan_cas9_prefix_hits_raw_range` wrapper remains the parity reference.

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

Finite early-stopping limits activate chunk-local guide counters. Workers mask
inactive guide bits before Myers verification and stop claiming chunks when all
guides are inactive. A guide retires only when a `(limit + 1)`th accepted hit
proves one exact-distance bucket incomplete. Detail output retains any valid
capped subset; its seven-column schema is unchanged.

#### Early-stopping benchmark result

GRCh38 benchmarks selected chunk-local reduction as the production design. For
61 Cas9 guides at distance 3 and 24 threads, prefixHash-style caps reduced
guide/window verification pairs from 1,583,279 to 122,692 and retired 51 guides;
detail median improved from 1.418 s to 1.359 s. An 11-run, single-thread Cas9
distance-4 count test measured a 2.27x paired speedup (84.94 s versus 36.80 s
median). Unlimited and default one-million caps remained near baseline. All
completed correctness comparisons passed.

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
    elseif geometry == Cas12a:
        simd_find_unambiguous_TTTV_and_BAAA_windows(bases)
    else:
        generated_simd_match(matcher_type, bases)
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

| Stage | Optimized canonical/generic d0-d4 paths | General `:legacy` path |
|---|---|---|
| Supported query | Eligible motifs, d0-d4, 16-base hash, <=64 guides per batch | Other motif sizes and hash lengths |
| Query structure | Presence bitmap + compact sorted directory + `UInt64` guide masks | `Dict` from hash to guide mask or guide-index vectors |
| Reference access | FAI range reads into reusable raw buffers | FASTA/2bit records loaded and converted to `LongDNA` |
| PAM search | Canonical or compile-time generic x86 SIMD masks evaluate 64 starts | `findguides` over materialized chromosome sequences |
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

### Generic-kernel cost relative to specialized Cas9

An August 2, 2026 scanner microbenchmark forced the canonical Cas9 motif through
the generic geometry so both kernels processed identical windows. It used 32 MB
of deterministic random A/C/G/T reference, distance 3, one thread, an empty
`UInt32` hash query, and 11 alternating timed runs after warmup. This isolates
motif scanning, prefix packing, and failed hash lookup; it excludes path/query
construction, verification, I/O, and result writing.

| Cas9 scanner | Median | Throughput | Relative latency |
|---|---:|---:|---:|
| Hand-written Cas9 | 66.7 ms | 480 MB/s | 1.000x |
| Forced typed generic | 78.6 ms | 407 MB/s | 1.178x |

Both kernels found exactly 4,001,477 motif candidates. The generic kernel took
17.8% longer and delivered 15.1% lower throughput, or about 85% of specialized
Cas9 throughput. This is a hot-loop result, not an end-to-end guarantee. Shared
query construction, verification, I/O, and output usually reduce the relative
effect; different PAM frequencies change candidate work.

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

1. The fastest kernels remain separate Cas9 and Cas12a geometries at distances
   0 through 4 with a 16-base prefix. Other eligible motifs use a typed generic
   SIMD geometry; unsupported sizes use the exact legacy path.
   Shared validation, bounds, overlap, and scheduling use
   `PrefixScanGeometry{Kind}` without moving motif branches into SIMD loops.
2. Each optimized query holds at most 64 guides. Larger public API and CLI
   searches rescan the reference once per sequential batch.
3. AVX2/BMI2 and AVX-512F/BW/BMI2 SIMD backends are available. Unsupported CPUs
   use the portable fused-directory or legacy path. There is no ARM kernel.
4. FASTA requires `.fai`; `.2bit` is streamed directly without a sidecar index.
5. Ambiguous query guides are rejected.
6. `ambig_max` supports zero through three IUPAC-ambiguous reference positions
   per complete guide/PAM window. Larger ambiguity allowances are unsupported.
7. Early stopping cannot cancel chunks already claimed by workers. Chunk-local
   reduction bounds this overshoot.
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
14. The exported three-argument `search_prefixHashScan` accepts registered names
    or custom Julia `Motif` objects at distances 0 through 4. The standalone CLI
    accepts registered names or complete custom motif definitions, including
    PAM position, strand subsets, and extension direction.
15. Source separates constant-free helpers, Cas9 and Cas12a kernels,
    verification, streaming, and orchestration in the parent `CHOPOFF` module.

## Potential speed optimizations

Unfinished items in this section are optional research directions. They do not
precede product completion, qualification maintenance, or portability work.

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
13. **Completed:** add an AVX-512F/BW profiling kernel that handles 64 input
    bytes per load and retains BMI2 prefix packing.

### Lower-priority micro-optimizations

14. Reduce FASTA newline-compaction copies or scan an mmap-backed representation
    directly. This adds format complexity and should follow I/O measurements.
15. Avoid string conversion and generic dedup keys during final commit. This is
    low priority for the 25,826-row Cas9 workload but credible for the
    364,581-row Cas12a workload.

## Feature status

Performance work should not be mixed blindly with generalization. Current
feature status is:

1. **Completed:** export and document registered and custom motifs at distances
   0 through 4 for `search_prefixHashScan`, add a standalone CLI search mode,
   and report the resolved backend and execution modes without enabling
   statistics.
2. **Completed:** support more than 64 guides through transparent sequential
   64-guide batching after the one-pass alternatives lost the GRCh38 benchmark.
3. **Completed:** computational early stopping masks retired guides and cancels
   future chunk claims. Detail subsets may vary with scheduling by design.
4. **Completed:** add distance 4 at p16 as the functional upper bound and
   benchmark mode.
5. **Completed:** add Cas12a as a separate specialized scan geometry using its
   precomputed paths without forcing Cas9 constants into a generic SIMD loop.
6. **Completed:** add optimized 2bit streaming and bounded IUPAC-reference
   behavior for `ambig_max=0:3`.
7. **Completed:** add count-only output without traceback or detail rows.
8. **Completed:** add complete custom motif definitions to the standalone CLI,
   including PAM-left, PAM-right, internal-PAM, PAMless, strand-subset, and
   extension-direction configurations.
9. **Completed:** run the representative generic qualification matrix against
   prefixHashDB. It covers full GRCh38 Cas9-NGA, CasX, and 25-base-guide cases,
   65-guide multi-batch searches, distances 0 through 4, ambiguity zero through
   three, and bounded internal-PAM, PAMless, 16-base-guide, strand-subset,
   FASTA/2bit, IUPAC, indel, and chunk-boundary cases. All 220 detail cases
   passed the reference-backed parity classifier and all 220 count comparisons
   passed. Of the detail cases, 147 were exact; the remainder contained only
   classified prefixHashDB ambiguity-limit or duplicate-row behavior.
10. **Completed:** make multi-batch detail output atomic. Batches append to a
    sibling temporary file; successful completion renames it over the requested
    path, while failure removes it and preserves any previous output.
11. **Completed:** add typed AVX-512F/BW profiling with explicit Julia and CLI
    selection, CPU-qualified automatic dispatch, parity tests, and ZMM codegen
    verification.

## Remaining implementation priorities

The following decisions are closed:

- Distance 4 remains fixed at p16. The p14/p15/p16 work selected p16, and d4 is
  intentionally a stretch configuration rather than an optimization target.
- Sequential 64-guide batching remains the large-guide architecture. It was
  the fastest of the tested multi-guide designs; no new one-pass representation
  or automatic routing policy is planned.
- Count-only output and direct prefixHashDB parity are complete.

Remaining work, in priority order:

1. **Product completion.** Add clearer path/query memory and progress reporting.
2. **Qualification maintenance.** Add randomized property tests and extend
   long-guide coverage beyond 28 bases using a suitable exact oracle. The
   completed prefixHashDB matrix cannot cover those d4 candidates because its
   packed representation is limited to 32 bases.
3. **Portable performance.** Qualify scalar and fused fallbacks on CPUs without
   AVX2/BMI2, then add an ARM SIMD path if profiling justifies it.

Further Cas9/Cas12a micro-optimization, NUMA experiments, alternative
prefilters, mmap-backed FASTA, and vectorized traceback remain
lower-priority research directions. Continue one only with exact parity and a
credible end-to-end gain.

## Recommended decision rule

Continue speed research while an experiment has a credible route to at least a
10% end-to-end improvement on GRCh38. No-statistics execution, scan/verify
fusion, buffered-vector reuse, global scheduling, parameter sweeps, 2 MiB
chunks, parallel guide hashing, production profiling, and radix-ordered lookup
are now measured.

Cache-local directory lookup passed, but its paired end-to-end gain was 6.2-6.7%
and lookup is only 7.9% of cumulative samples at 24 CPUs. Cas12a specialization
is now complete and confirms that motif geometries should remain separate.
If performance research resumes, Cas9 work should target genome-order
presence-bitmap probes or explain worker wait through NUMA/pinning measurements.
Cas12a detail work should target accepted-candidate verification, alignment
materialization, deduplication, and commit. Require exact parity and report both
motifs so an optimization for the result-heavy Cas12a workload does not regress
Cas9.

The likely remaining gain without a genome index is meaningful but smaller than
the previous 10x improvement. Another 1.5-3x is plausible only if profiling
confirms avoidable lookup, scheduling, or temporary-data costs. Another 10x is
unlikely because every exact indexless search must still inspect the reference.

## Open research questions

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
- one bounded reference scan per 64-guide batch;
- deterministic results with prefixHashDB-compatible coordinates and distances;
- specialized hot loops for important motif geometries, without motif branches
  inside those loops.

The useful output surface, computational early stopping, custom CLI motif
definitions, and representative generic correctness qualification are now
complete. Atomic multi-batch detail output is also complete. Operational
reporting is the next implementation priority.

Implementation order changed July 22, 2026: detailed output for distances 0
through 3 shipped first at p16, followed by distance 4 as the functional and
benchmarking ceiling.

### Priority 1: production distances 0 through 4 (completed at p16)

Distances 0, 1, 2, 3, and 4 are supported public configurations for registered
and custom motifs. A lower requested distance does not run a larger-threshold
query and discard higher-distance results. Guide lengths 20 and 21 reuse exact
p16 assets; other lengths generate the requested paths once per search.

Distance parameterizes path selection, chunk overlap, and the Myers threshold.
The existing Cas9 and Cas12a PAM/profile kernels remain separate and are not
copied once per distance. Dispatch resolves the distance before entering each
hot loop. Production prefix length remains fixed at 16.

The p14/p15/p16 evaluation retained p16 as the final production prefix length.
Distance 4 is supported as a stretch configuration, not as an optimization
target. Exact detail and count parity with prefixHashDB remain release gates at
every supported distance.

### Priority 2: count output mode (completed)

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

Finite limits use chunk-local counters. A guide retires after a `(cap + 1)`th
hit proves incompleteness in any exact-distance bucket. Subsequent candidates
mask that guide before Myers verification, and workers stop claiming chunks
when every guide retires. Count output caps each bucket and reports
`complete=false`; non-triggering buckets may then be partial lower bounds.

This mode is also a performance feature. In the current 61-guide human
experiments, detail output performs 25,826 Cas9 tracebacks and 364,581 Cas12a
tracebacks. Count mode bypasses both workloads while preserving exact
per-distance counts.

An August 8, 2026 GRCh38 check used 61 guides, 24 Julia threads, one warmup,
and three alternating timed runs per mode. Count/detail summaries had exact
parity and count diagnostics performed zero tracebacks:

| Motif | Distance | Detail median | Count median | Speedup |
|---|---:|---:|---:|---:|
| Cas9 | 3 | 1.556 s | 1.446 s | 1.08x |
| Cas9 | 4 | 14.376 s | 12.556 s | 1.14x |
| Cas12a | 3 | 3.819 s | 1.069 s | 3.57x |
| Cas12a | 4 | 20.819 s | 10.758 s | 1.94x |

### Completed: thousands of guides through batching

The `UInt64` guide mask remains the fast representation for every optimized
query. Larger input sets use sequential 64-guide batches, which bounds live
query and detail-result state while preserving exact per-guide behavior.

The rejected one-pass experiment used a compact hash-to-guide-ID directory and
optional reference-aware query filtering. Although filtering reduced the
dispersed query from roughly 1.45 GB to 242 MB, full-detail latency remained
about three times slower than batching. The dominant workload contains tens of
millions of verified detail rows, so retaining one monolithic result and dedup
state outweighed the saved reference scans.

Sequential batching is the production architecture rather than a fallback.
The tested one-pass and large-directory alternatives were slower, so large-guide
query redesign is not active work.

### Completed: bounded IUPAC ambiguity

`motif.ambig_max=0:3` bounds ambiguous reference positions across the complete
candidate window. PAM symbols use IUPAC compatibility masks. Ambiguous hashed
prefixes expand to compatible concrete hashes, and ambiguity outside the prefix
uses IUPAC-aware Myers verification and traceback. Query guides remain
unambiguous. FASTA preserves supported IUPAC symbols; 2bit represents ambiguous
blocks as `N`.

### Completed: distance 4 benchmark ceiling

Distance 4 is the largest public search threshold and a benchmarking ceiling.
For 20- and 21-base guides at p16, it uses the existing split d4 matrices:
8,196,801 symbolic paths and about 125 MiB. Other guide lengths generate and
deduplicate their paths before scanning.

This mode intentionally reuses the d0-d4 compact query and motif-specific scan
architecture. It is not expected to match lower-distance speed or memory use. In
an earlier isolated 61-guide Cas9 memory run with 24 query workers, query
construction took 10.8 seconds, produced 111,720,240 guide/hash associations
and a 1.09 GB final query, reaching 4.51 GB peak process RSS. No guide-count cap
below the existing 64 is imposed. The p14/p15/p16 evaluation retained p16;
compressed or staged d4 representations are out of scope because d4 is already
at the practical edge of the algorithm.

### Workload-specific performance after generalization

After the core implementation priorities, Cas9 and Cas12a performance research
should retain separate gates:

- Cas9 work should investigate genome-order presence-bitmap latency,
  last-level-cache behavior, NUMA-local query copies, and scheduler tail
  imbalance.
- Cas12a detail work should target accepted-hit materialization, traceback,
  string construction, deduplication, and CSV commit.
- Count mode should be benchmarked independently because it removes most of the
  Cas12a-specific output cost and exposes the scan/verification ceiling.
- Non-AVX2 fallback qualification and ARM SIMD belong to the portability
  priority. Mmap-backed input and further prefilter experiments remain optional
  speed research.

Continue an optimization only when it has a credible route to at least 10%
end-to-end improvement in its intended workload. Always report detail and
count modes separately and verify that a motif-specific win does not regress
the other production geometry.

### Successor readiness relative to prefixHashDB

For a one-shot, complete search with at most 64 guides on supported x86
hardware, `prefixHashScan` is already the likely successor to `prefixHashDB`.
It avoids database construction and storage, and the measured Cas9/Cas12a
searches are substantially faster. It cannot yet be described as universally
superior because the two algorithms amortize genome work differently.

`prefixHashScan` scans the reference once for each 64-guide batch. Thus 61,
1,024, and 4,096 guides require 1, 16, and 64 reference scans respectively.
`prefixHashDB` pays genome processing once during database construction and can
reuse that index for arbitrary future searches. A persistent index can therefore
win for repeated searches even when one scan is faster than one indexed search.
The architectural comparison remains:

```text
prefixHashDB total = database build + searches * indexed search
prefixHashScan total = searches * ceil(guides / 64) * reference scan
```

The main remaining gaps are:

1. **Product completion.** Add path/query memory and progress reporting.
   Custom motif definitions in the CLI and atomic multi-batch detail output are
   complete.
2. **Qualification maintenance.** The representative generic matrix is
   complete: 220 detail cases passed its reference-backed classifier and all
   220 count comparisons passed. Add randomized property tests and qualify
   guide lengths above 28 with an oracle other than prefixHashDB, whose packed
   d4 representation cannot cover those candidates.
3. **Portable performance.** Qualified AVX-512 and AVX2/BMI2 systems use the
   fastest streaming path.
   Correct fused-directory and legacy fallbacks exist, but must be qualified
   against prefixHashDB on unsupported CPUs. ARM SIMD is the primary missing
   optimized backend.

The generic kernel's measured 17.8% Cas9 hot-loop latency penalty is not a
replacement blocker. The decisive remaining product limitation is incomplete
operational reporting. Sequential batching, atomic multi-batch output,
computational early stopping, d4/p16, custom CLI motifs, and representative
generic qualification are closed implementation items.

#### Replacement gates

Make `prefixHashScan` the documented default for its qualified workload when:

1. Canonical searches retain exact detail parity, and representative generic
   distances 0 through 4 retain exact or reference-classified parity on sample
   and human-scale fixtures. Randomized coverage and guide lengths above 28 are
   tracked as qualification maintenance.
2. Scalar, fused, FASTA SIMD, and 2bit SIMD backends have identical results;
   unsupported configurations select a correct fallback rather than fail.
3. Count output marks early-stopped rows incomplete; detail output documents its
   scheduling-dependent valid subset. Cancellation reduces verification or
   unclaimed work.
4. Progress/memory reporting is clear enough for production use. Atomic
   multi-batch detail output and custom CLI motifs are already supported.
5. Non-AVX2 fallbacks are qualified, with ARM SIMD tracked as the primary
   portability extension.
6. Cas9/d3 and Cas12a/d3 detail latency regresses by no more than 3% unless a
   measured feature-level benefit justifies it.

The initial default should be one-shot searches with eligible motifs. Larger
guide sets continue through the measured sequential batching path. Keep
`prefixHashDB` as the persistent-index backend for repeated or heavily capped
workloads; no automatic crossover policy is currently planned.
