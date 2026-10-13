# Variant-aware prefixHashScan

`search_prefixHashScan` accepts a normalized VCF through `vcf` or `--vcf`.
The first release searches sequence-resolved SNVs, MNPs, insertions, deletions,
and unequal-length replacements. It searches compatible REF/ALT combinations up to
an ALT-allele cap. It does not interpret genotypes, phase, FILTER, or allele frequency.

Status: complete for single-sample VCFs, validated on GIAB HG001 genome-wide.
Cohort and catalog VCFs (gnomAD, ClinVar) work per chromosome with a small cap;
genome-wide cohort runs need the deferred memory work below.

A search covers the non-conflicting allele combinations represented by accepted
records, with at most `vcf_max_alleles` (default 3) ALT alleles in one local
haplotype, subject to the motif's reference-ambiguity policy. This is a union of
possible sequences. It does not establish that a sample carries them. Finite
early-stopping limits return a subset. Exhaustive combination is infeasible for
cohort and catalog VCFs; see [Allele cap and measured cost](#allele-cap-and-measured-cost).

## Implementation and validation checklist

Baseline: the original design in revision `4523238e`, amended by the agreed v1
scope. Implementation and validation are separate: a checked implementation can
still lack targeted tests. “Tested” means the named bounded cases passed, not that
every possible genome has been exhausted. Deferred features are not implemented
merely because their inputs have rejection tests.

| ID | Requirement | Implementation | Validation / evidence |
|---|---|---|---|
| V01 | Sorted records, contig names, bounds and REF validation | [x] implemented, [reader](src/variants/vcf_read.jl) | tested; [isolated input failures](test/variants/input.jl) |
| V02 | Plain VCF, gzip and BGZF | [x] implemented, reader | tested; [byte parity, genuine BGZF split within a row](test/variants/input.jl) |
| V03 | Ignored markers and unsupported-allele errors | [x] implemented, reader | tested; [marker siblings, aggregate warning, rejection fixtures](test/variants/input.jl) |
| V04 | SNV/MNP/indel/replacement combinations and reference choices | [x] implemented, [enumerator](src/variants/enumerate.jl) | tested; [full-haplotype oracle and required witnesses](test/variants/corpus.jl) |
| V05 | Pairwise conflicts and original-record ALT exclusivity | [x] implemented, reader/enumerator | tested; [same-gap, boundary, nontransitive and split-source cases](test/sample_data/vcf/README.md) |
| V06 | IUPAC ALT expansion (at most 64 realizations per ALT) and reference ambiguity policy | [x] implemented, reader/scanner | tested; [all symbols](test/variants/corpus.jl), [ambiguity budgets 0–3, packed fragments, expansion limit](test/variants/coverage.jl), [REF compatibility](test/variants/input.jl) |
| V07 | Window reach, deletion bridges, chromosome/chunk edges | [x] implemented, enumerator/[scanner](src/variants/search.jl) | tested; [spacing 1 through W+1](test/variants/boundaries.jl), [bridge and edge fixtures](test/sample_data/vcf/README.md) |
| V08 | PAM changes, strands, motifs and boundary bulges | [x] implemented, scanner | tested; [creation/destruction and generic motifs](test/variants/corpus.jl), [Cas9/Cas12a, both strands, distances 0–4](test/variants/boundaries.jl) |
| V09 | Endpoint lifting and insertion-origin axes | [x] implemented, [fragments](src/variants/fragments.jl) | tested; [full-row oracle and inserted endpoints](test/variants/corpus.jl), [long-fragment execution matrix](test/variants/execution.jl) |
| V10 | Complete allele and INFO/FORMAT annotations | [x] implemented, reader/fragments | tested; [A/R/G/fixed/variable cardinalities, multiple IDs, escaping and missing values](test/variants/input.jl) |
| V11 | Reference-overlap annotations | [x] implemented, fragments | tested; [complete-row oracle](test/helpers/vcf.jl), [exact multiple-ID annotation](test/variants/input.jl) |
| V12 | Independent caps; duplicates and flank candidates spend no budget | [x] implemented, scanner | tested; [cap retirement and duplicate-fragment regressions](test/variants/regressions.jl) |
| V13 | Preserve allele combinations in filtering; summarize loci once; insertion starts within `distance` of an insertion end join POS | [x] implemented, [utilities](src/utils.jl) | tested; [exact summaries, origin axes, ordering and idempotence](test/variants/execution.jl) |
| V14 | Region selection, API/CLI parity and atomic output | [x] implemented, [API](src/db_prefix_hash_scan.jl) | tested; [input failures](test/variants/input.jl), [CLI parity and a failure injected into the second batch's variant scan](test/variants/execution.jl), [shared atomic rollback helper](test/src/prefix_hash_scan.jl) |
| V15 | Batching, reference formats, SIMD and thread parity | [x] implemented, API/scanner | tested; [64/65/129 guides, FASTA/2bit, chunks, available SIMD and workers](test/variants/execution.jl) |
| V16 | Reproducible permanent complex fixtures | [x] implemented | tested; [42 sequence cases, 2 metadata/marker cases, 14 invalid inputs](test/sample_data/vcf/README.md), byte-for-byte regeneration |
| V17 | Seeded/dense extended suite with reproducible failures | [x] implemented | tested; [200 seeded cases, 60 dense mixed cases and 6/8/10-variant exhaustive combinations](test/extended/variants.jl) |
| V18 | Sparse chromosome-scale timing | [x] implemented, [benchmark](scripts/benchmark_prefix_hash_scan_vcf.jl) | tested on synthetic chr22; measurements below |
| V19 | `vcf_max_alleles` cap, truncation and above-default warnings | [x] implemented, enumerator | tested; [capped oracle, warnings, validation](test/variants/coverage.jl), [capped dense cases](test/extended/variants.jl); chr22 cost table below |
| V20 | Realistic inputs: GIAB HG001, ClinVar, gnomAD v3.1.2 chr22 | [x] validated | windowed oracle, reference parity, byte parity and an independent row validator; see [Realistic validation](#realistic-validation) |
| D01 | Genotypes, phase blocks and carrier/sample summaries | [ ] deferred | untested; v1 ignores GT/phase |
| D02 | Symbolic SVs, BNDs and nearby small-variant combinations | [ ] deferred | untested; v1 rejects non-literal alleles |
| D03 | Interval/BED selection with deletion-aware reachability | [ ] deferred | untested; whole chromosomes only |
| D04 | IUPAC compression valve | [ ] deferred | untested; explicit expansion only |
| D05 | `lost` summaries | [ ] deferred | untested; absence of a row is not proof of loss |
| D06 | Legacy vcfDB retirement | [ ] deferred | legacy implementation/tests remain |
| D07 | Parallel, per-chromosome preparation | [ ] deferred | preparation is serial and holds all selected chromosomes; run cohort VCFs one chromosome at a time with `region` |
| D08 | Representative population and patient-scale measurements | [~] partial | chr22 cost table for GIAB, ClinVar and gnomAD tiers; genome-wide cohort timing untested |
| D09 | Variant-phase speed (see [Deferred design](#deferred-design)) | [~] partial | work sizing, stored reference spans and block reads done; see [Speed changes](#speed-changes) |

Deliberate changes from the original design: unknown contigs and REF mismatches
fail immediately; filtering preserves each allele combination; inserted endpoints
retain an origin axis and offset; reference and variant work use separate scheduling
phases. The full original roadmap is not a claim about delivered v1 behavior.

## Run a search

Normalize and sort the VCF against the exact search reference. Contig names must
match the reference (`chr22` versus `22`); rename with `bcftools annotate
--rename-chrs` and drop contigs absent from the reference. Joining multiallelic
records (`-m+any`) keeps their ALT choices exclusive:

```bash
bcftools norm -f reference.fa -c e -m+any input.vcf.gz -Oz -o normalized.vcf.gz
bcftools sort normalized.vcf.gz -Oz -o sorted.vcf.gz

CHOPOFF search --guides guides.txt --output detail.csv --distance 3 prefixHashScan \
    --genome reference.fa --vcf sorted.vcf.gz --region chr22 --vcf_max_alleles 3
CHOPOFF filter --detail_file detail.csv --output filtered.csv --distance 3
CHOPOFF summarize --detail_file filtered.csv --output counts.csv --variant_aware --distance 3
```

Keep multiallelic records together when possible. The reader splits them internally
and remembers that their ALT choices are exclusive. Pre-split records without
source provenance become independent choices, subject to footprint conflicts.
That can add combinations absent from the original multiallelic record.

A tabix index is unnecessary. The reader streams plain VCF or gzipped VCF,
including concatenated BGZF members, through the existing decompressor.
BCF input is outside this release.

Use upstream tools for FILTER, frequency, and sample selection. Selecting sample
columns alone does not remove unused cohort ALT alleles. CHOPOFF searches every
retained ALT regardless of GT. See the
[bcftools filtering and allele-trimming options](https://samtools.github.io/bcftools/bcftools.html).

The Julia API extends both existing entrypoints:

```julia
search_prefixHashScan(guides, genome_path, output_file;
    motif = "Cas9", distance = 3,
    vcf = "sorted.vcf.gz",
    vcf_annotate = ["INFO/AF"],
    vcf_max_alleles = 3,
    region = "chr1,chr22")

search_prefixHashScan(guides, genome_path, motif, output_file;
    vcf = "sorted.vcf.gz", distance = 3, stats = PrefixHashScanStats())
```

`region=nothing` searches every reference chromosome. Otherwise, `region` accepts
comma-separated chromosome names. Unknown names, interval strings, and BED paths
are errors. Chromosome selection also works without a VCF, including count output
and the legacy scanner.

VCF searches require detail output and an eligible streaming geometry.
Invalid configuration fails before output opens. The public three-argument method
prepares variants once, reuses them across guide batches, and replaces output
atomically after all batches succeed.

## Input contract

The reader checks coordinate order, contiguous chromosome blocks, chromosome names,
reference bounds, and REF sequence. A mismatch fails the run. REF ambiguity must
contain every base allowed by the supplied reference at that position.
It does not create alternative reference backbones.

For chromosome-restricted runs, the reader checks ordering and chromosome names
throughout the file. It checks REF and ALT content on selected chromosomes.

ALT accepts DNA IUPAC symbols in either case. `R` becomes separate A and G choices.
`N` becomes A, C, G, and T. Each output alignment uses concrete ALT bases and retains
the original allele annotation. An ALT with more than 64 concrete realizations
(for example `NNNN`) encodes unknown sequence; the reader ignores it with an
aggregate warning. ClinVar contains such ALTs with up to 29 `N`s. Reference
ambiguity in flanks still uses `ambig_max`, which defaults to zero.

`ALT=.` has no alternate sequence. The reader also ignores spanning-deletion `*`
and gVCF `<NON_REF>`/`<*>` markers, with aggregate warning counts. Other ALT alleles
on the same record remain available. These markers do not supply a deletion or an
unknown alternate sequence by themselves.

Contig and REF errors name a bcftools fix: `annotate --rename-chrs` when the
reference has the same name with or without `chr`, `view -t ^contig` for contigs
absent from the reference, and `norm -f` for REF mismatches.

Symbolic SVs and breakends fail with an unsupported-allele error. Literal sequence
alleles have no size cutoff. Users must left-align upstream. The reader trims common
padding but does not validate left alignment or merge equivalent records.

Use `vcf_annotate` or `--vcf_annotate INFO/AF,FORMAT/AF` to retain selected fields.
Unknown header fields fail at startup. INFO/FORMAT fields with Number=A select the
chosen ALT value. Number=R retains the REF and chosen ALT values. Other cardinalities
retain their original value. FORMAT annotations name their samples without interpreting GT.

## Walks, fragments, and coordinates

A replacement occupies a half-open reference interval `[first, stop)`.
An insertion occupies the gap before `first`, with `first == stop`.
Base intervals conflict if they intersect. An insertion conflicts with an interval
only when its gap lies strictly inside that interval. Insertions at the same gap
conflict. ALT choices from one source record also conflict.

Conflict is pairwise. A deletion overlapping two mutually compatible variants does
not make those two variants exclusive. A SNP and an insertion immediately after
it can coexist.

Each concrete ALT anchors a local depth-first walk. The walk includes a reference
left flank of at most `W = motif span + distance` bases. It branches over subsequent
compatible variants until it emits W bases beyond the anchor. Deletions can therefore
bring distant reference positions into the same site. Skipping a record leaves the
sequence cursor unchanged. Insertions precede replacements starting at the same gap.
A walk takes at most `vcf_max_alleles` ALT alleles, counting the anchor. When a
compatible variant in reach is skipped only because of the cap, the anchor counts
as truncated.

Fragments use flat segment arrays. A segment either refers to a reference range or
to a concrete ALT. Empty ALT segments retain deletion junctions. Descriptors persist
across guide batches. Preparation reads the merged reference spans that fragments
use once, through a forward 1 MiB block reader, and keeps them in memory. Workers
materialize fragments from these spans into reusable buffers without file reads,
then use the existing prefix kernel, Myers verifier, and traceback.

Variant work packs fragments into items of about a quarter of the fragment bases
per scan thread (at least 64 KiB, at most the chunk size, 2 MiB by default), with
`distance + 1` gap bytes between fragments. Large fragments split into overlapping scan chunks.
Candidate windows must belong to one fragment and overlap its anchor before Myers.
After traceback, the consumed target plus PAM determines the allele annotations.
Unused extension bases do not contribute annotations.

A best alignment can avoid the anchor because nearby variants change candidate
context. Keep that alignment when it touches another taken allele. Exact-record
deduplication removes repetitions from other fragments before counting them.
Rejecting such alignments by their leftmost consumed allele loses valid rows at
chromosome edges.

Reference work and variant work run as separate scheduling phases with independent
budgets. They share the query directory and guide profiles. Unique accepted variant
records consume the variant budget after site annotation and deduplication checks.
A retired reference guide does not retire its variant search. Each phase uses the
supplied `early_stopping` vector. Variant retirement produces an incomplete-output warning.
Uncapped results are deterministic across worker schedules.

VCF output retains the original seven columns and adds:

| Column | Meaning |
|---|---|
| `end` | Opposite endpoint of the consumed site, including PAM |
| `variants` | Applied alleles touching the site, separated by semicolons |
| `overlapping_variants` | Reference rows only: input alleles intersecting the site |
| `start_origin`, `end_origin` | ALT identity and 1-based trimmed ALT offset when the endpoint lacks a reference coordinate |

`start` retains the existing prefixHashScan convention. Substitution bases map to
corresponding reference positions. Inserted bases and unequal-length replacements
project to their VCF POS and retain an allele-relative origin. Their origin format
is `chrom:POS:REF>original_ALT:offset`. Concrete IUPAC realizations share that origin
axis. Distinct positions inside a long insertion remain distinct positions.

An allele annotation has the form `ID|chrom:POS:REF>original_ALT`. Missing IDs use
the allele coordinates. IUPAC expansions append `|resolved=concrete_ALT`.
Requested fields follow in request order. Tokens escape percent, semicolon, and pipe
as `%25`, `%3B`, and `%7C`. The CSV writer quotes fields containing commas or quotes.

A reference-equivalent ALT choice remains an annotated variant row. Absence of a
variant row does not prove that a variant destroys a site. No `lost` count or
cutting-activity prediction exists in this release.

## Filtering and summaries

Filtering assigns loci across reference and variant rows together. It uses the
existing distance-radius rule within each guide, chromosome, strand, and coordinate
axis. A start on an inserted base within `distance` of either insertion end joins
the reference axis at POS, so a short insertion does not split one locus in two.
Deeper inserted positions use their ALT origin axis and offset, so distinct sites
inside a long insertion stay distinct.

Within each locus, `filter_overlapping` retains the best row for every distinct
allele annotation, including the reference choice. It preserves different allele
combinations even when their sequences match. On a distance tie, a reference row
represents the locus for summary counts.

`summarize_offtargets(...; variant_aware=true, distance=k)` counts each locus once
at its best distance. It accepts raw or filtered VCF output. Use the same explicit
`distance` for filtering and summarizing. Plain summaries retain their row-count
behavior. Searches without `vcf` retain the original output schema.

## Validation and measurement

The regular entrypoint is [test/src/variants.jl](test/src/variants.jl). Its
[oracle](test/helpers/vcf.jl) constructs complete alternative haplotypes with a
separate edit model, searches them with the brute-force legacy backend, and lifts
coordinates with its own per-base map. It compares all twelve output columns,
including allele annotations, reference overlaps and insertion origins, and rejects
duplicate rows. Hand-calculated checks validate the oracle's edit model. Production
VCF parsing, enumeration, conflict checks and coordinate lifting do not participate
in the oracle; alignment and traceback remain shared with the legacy search.
[test/src/distance_metrics.jl](test/src/distance_metrics.jl) checks `align` against
a naive edit-distance DP on 5,000 seeded guide/reference pairs with k = 1-4,
including IUPAC reference bases and indels; distances, gapless sequences and
column costs all agree.

[Permanent fixtures](test/sample_data/vcf/README.md) combine the existing
semirandom genome example with tiny companion references. Every positive sequence
case requires a specific allele witness, so an empty result cannot pass. Invalid
inputs are separate files. The generator's `--check` mode regenerates into a
temporary directory and checks bytes and file lists without changing the corpus.

The regular suite includes exhaustive bounded haplotypes, motif/strand/distance
sweeps and the public API/CLI contract. The separate extended suite uses a fixed
seed, up to eight source records per random case, and dense combinations up to
1,024 haplotypes. On row mismatches, both suites preserve inputs and differences
outside temporary-directory cleanup. Extended failures also retain the seed and
reproduction command. These bounded tests do not establish population-scale cost
or validate deferred features.

From the package root:

```bash
source scripts/dev_env.sh
julia --project=. -t 4 -e 'include("test/src/variants.jl")'
julia --project=. -t 4 test/extended/variants.jl
julia --project=. test/sample_data/vcf/make_test_vcfs.jl --check
(cd test && julia --project=.. -t 4 runtests.jl)
# Optional alternative seed or larger random sample:
CHOPOFF_VCF_SEED=123 CHOPOFF_VCF_CASES=1000 julia --project=. -t 4 test/extended/variants.jl
```

[Coverage tests](test/variants/coverage.jl) add multiple chromosomes (including
an empty middle chromosome and an edge fragment packed after another), several
guides on random backgrounds, ambiguity budgets with packed fragments, REF/ALT
trimming paths (shared suffix, multiallelic STR, interior MNP, IUPAC insertion),
insertion gaps at site boundaries, Cas9/Cas12a edges on both strands, soft-masked
references, per-guide variant early stopping, the IUPAC expansion limit and the
ALT cap against a capped oracle (haplotypes with at most k edits). Runtime
mutations of the overlap rule, coordinate lifting, suffix trimming and the cap
each make these tests fail. The extended suite adds 60 dense mixed cases with
indels, MNPs, IUPAC and multiallelic records, every fourth with a random cap.

Local validation on 2026-10-12 used Julia 1.10.10 with four threads. The regular
suite passed 3,492 assertions, including fixture regeneration; the extended suite
passed 1,195 with seed `20261009`. The full package suite also passed with four
threads. The execution matrix exercised portable, auto,
AVX2 and AVX-512 backends, one/four workers, FASTA wrapping at 7/60 bases, 2bit,
and chunks of 64, 97 and 2,097,152 bases. These are correctness checks, not timings
for a representative biological dataset.

The annotation regression exposed cumulative ID ranges in the installed VCF parser.
The reader preserves the complete ID field so multiple IDs appear exactly once in
both variant and reference-overlap annotations.

`PrefixHashScanStats` records input counts, fragment counts and bases, the largest
fragment count of one anchor (`max_anchor_fragments`), anchors truncated by the cap
(`vcf_truncated_anchors`), descriptor bytes (including stored reference spans),
preparation time, variant fragment materialization time (`variant_io_ns`), candidates, anchor rejections, accepted variant rows, and
retired variant guides. `verbose=true` reports preparation
counts and descriptor memory once before guide batching. It also logs throttled
`prefixHashScan progress` records (every `progress_interval` seconds) with
`batch` set to `VCF read`, `VCF enumerate` or `<batch> VCF` for the variant scan.
The CLI always runs verbose.

Report performance as measured wall time and memory. Additional scanned bases alone
do not predict overhead: small reference reads, enumeration, traceback, and output
can dominate. Real-VCF costs are in [Allele cap and measured cost](#allele-cap-and-measured-cost).

Run `scripts/benchmark_prefix_hash_scan_vcf.jl GENOME VCF GUIDES CHROM|all OUTDIR [RUNS]`
with `julia --project=. -t 8`. It warms both modes, alternates measured runs,
checks reference-row parity and repeatability, and writes timings and allocation
counts. Peak RSS includes warmup and both modes. It does not measure per-mode
live memory.

A local GRCh38 chromosome-22 check used 61 Cas9 guides, distance 3, eight threads,
and five alternating runs. The synthetic dataset contained 62,686 records at least
600 bases apart: eight SNPs, one single-base insertion, and one single-base
deletion per ten records, skipping ambiguous reference positions. The medians were
0.457 s for reference alone and 1.609 s with VCF (3.52 times reference alone).
The VCF scan built 62,686 fragments covering 3.32 million bases, with 25.1 MB of
descriptors. Preparation took a median 0.524 s. Worker I/O totaled 0.366 s.
The runs allocated about 310 MB and 643 MB, respectively. Process peak RSS reached
1.29 GB. All reference rows matched, and repeated outputs were byte-identical.
These sparse synthetic records do not measure dense population combinations.
The figures above predate the speed changes below.

### Speed changes

Deferred-design items 1-3 were measured on October 13, 2026 with the same
benchmark: GRCh38 chromosome 22, 61 Cas9 guides, distance 3, eight threads, five
alternating runs, default `vcf_max_alleles` 3. Worker time is summed over workers.
Variant output was byte-identical to the previous code for all three inputs.

| VCF | Before | After | Preparation before / after | Worker materialization before / after |
|---|---:|---:|---:|---:|
| Synthetic (62,686 records) | 1.56 s | 1.15 s | 0.49 / 0.36 s | 0.38 / 0.03 s |
| GIAB HG001 v4.2.1 | 8.10 s | 7.38 s | 6.92 / 6.56 s | 0.36 / 0.03 s |
| ClinVar, multiallelics joined | 74.5 s | 32.8 s | 30.7 / 30.3 s | 352 / 7.6 s |

Reference-only medians were 0.44-0.49 s. Stored reference spans add 5 MB
(synthetic), 4 MB (GIAB) and 1 MB (ClinVar) of descriptors. Items per thread of
2, 4 and 8 gave the same synthetic median within noise; the code uses 4.

GIAB and ClinVar preparation was dominated by reading the whole genome-wide VCF,
since `region` still reads every record. Reading all chromosomes with REF checks
took 67.0 s before and 38.7 s after the block reader for GIAB, and 81.4 s and
52.4 s for ClinVar; reading with no chromosome selected took 9-13 s.

The reader now splits records by line and parses only CHROM and POS for
unselected chromosomes. Name, contiguity and order checks still cover every
record; the full parser runs on selected records only. Most of the old cost was
not the parser: the reader loop ran inside an `open(...) do` closure, which boxed
the reassigned loop state (`previous_chrom`, `previous_pos`). Without the closure,
reading GIAB with no chromosome selected takes 1.6 s instead of 6.3 s (line
splitting and decompression take 0.9 s). Medians of five alternating runs, same
benchmark, measured October 13, 2026, output byte-identical:

| VCF, chromosome 22 | Before | After | Preparation before / after |
|---|---:|---:|---:|
| Synthetic | 1.00 s | 0.90 s | 0.39 / 0.35 s |
| GIAB HG001 v4.2.1 | 6.90 s | 2.68 s | 6.41 / 2.15 s |
| ClinVar, multiallelics joined | 33.6 s | 28.0 s | 31.3 / 26.4 s |

Genome-wide reads are unchanged (43.6 s before, 44.3 s after for GIAB), since every
record is parsed either way.

Descriptor memory was measured with `Base.summarysize`, which tracks object
identity. On genome-wide GIAB it took 189-223 s, four times the preparation, in
every CLI run (the CLI logs it) and every run with `stats`. A linear count of
array and string bytes replaced it (0.3 s). It omits per-object headers, so it
reports about 20% less (1.35 GB versus 1.71 GB for genome-wide GIAB).

### Genome-wide preparation profile

Preparation of every chromosome, k = 3, Cas9 distance 3, three repetitions in one
process after a chr21 warm-up, measured October 13, 2026. "Read" is
`read_scan_variants` alone; "preparation" repeats the read and adds enumeration
and reference spans.

| VCF | Read | Read GC | Preparation | Preparation GC | Allocated | Fragments | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| GIAB HG001 v4.2.1 | 29-37 s | 10-18 s | 41-44 s | 15-17 s | 19.6 GB | 4.42M | 4.1 GB |
| ClinVar, multiallelics joined | 47-95 s | 14-43 s | 1,034-1,386 s | 114-145 s | 325 GB | 455M | 158 GB |

For GIAB, the read is about 80% of preparation. It allocates 14.3 GB, about 3.6 KB
per record. A CPU profile of one GIAB preparation (5,157 samples at 10 ms, GC
time counted at the allocating line):

| Part | Share |
|---|---:|
| `foreach_concrete_variant`: copies and IUPAC checks run for plain ACGT ALTs too | 21% |
| ID and allele-key strings: `"$chrom:$pos:$ref>$alt"` built twice, `variant_token` runs three `replace` calls | 20% |
| `VCF.Record` parse of selected records | 17% |
| Enumeration walk and reference spans | 17% |
| REF and ALT accessors, `uppercase`, REF check | 12% |
| Variant sort (tuple key ending in the annotation string) | 5% |

Genome-wide ClinVar at k = 3 is not practical: 455M fragments, 91 GB of descriptors
and 158 GB peak RSS, with enumeration taking 95% of the time. Catalog VCFs need
`region` per chromosome, as the cost table shows.

Repeated GIAB preparations in one process took 40.6-44.0 s, so the 84 s and 100 s
benchmark runs were not reproduced with preparation alone. Their cause (host load,
or GC state left by the previous search) is unknown.

Candidate fixes for the read, none applied yet:

1. Pass plain ACGT ALTs to the callback directly in `foreach_concrete_variant`.
2. Build the `chrom:pos:ref>alt` key once per ALT, and skip `replace` in
   `variant_token` when the string has no `%`, `;` or `|`.
3. Take ID, REF and ALT from the tab-split line and keep `VCF.Record` only for
   `vcf_annotate`.

Items 1 and 2 are small changes. They touch 41% of the profile, so an estimated
20-30% of GIAB genome-wide preparation, not yet measured. Item 3 adds a second field parser.

## Allele cap and measured cost

`vcf_max_alleles` (CLI `--vcf_max_alleles`, default 3) limits the ALT alleles in one
local haplotype. Every site whose local haplotype (window plus extension) carries at
most k ALT alleles is reported; combinations with more are not. The run warns once
with the number of truncated anchors and example loci. Setting k above the default
warns that memory can grow quickly.

Fragments per anchor grow like `sum(C(n, j) for j < k)`, where n is the number of
compatible variants in reach. Cohort and catalog VCFs place dozens of overlapping
indel alleles at STR and VNTR loci, so cost rises steeply with k. Single-sample VCFs
stay small at any k.

GRCh38 chromosome 22, 61 Cas9 guides, distance 3, eight threads, October 11, 2026
(reference-only search 3.5 s). Peak RSS is the process peak. Counts marked "count"
come from enumeration only, without a search.

| VCF (chr22 records) | k | Fragments | Peak RSS | Time | Truncated anchors |
|---|---:|---:|---:|---:|---:|
| GIAB HG001 v4.2.1 (49,484) | 2 | 55,680 | 0.96 GB | 21 s | 662 |
| | 3 | 56,994 | 0.93 GB | 21 s | 172 |
| | 4 | 57,580 | 0.93 GB | 22 s | 58 |
| | 8 | 58,047 (count) | | | 1 |
| ClinVar, multiallelics joined (86,900) | 2 | 1.0M | 1.2 GB | 27 s | 73,084 |
| | 3 | 12.5M | 4.5 GB | 91 s | 62,672 |
| | 4 | 155M (count) | about 50 GB | | 54,601 |
| | 5 | 1.6B (count) | | | 47,980 |
| gnomAD v3.1.2 PASS, AF ≥ 1% (231,346) | 2 | 3.4M | 2.9 GB | 38 s | 88,522 |
| | 3 | 30.4M | 10.1 GB | 187 s | 56,089 |
| | 4 | 252M (count) | about 80 GB | | 39,530 |
| gnomAD PASS, AF ≥ 0.1% (472,415) | 2 | 6.7M | 3.5 GB | 70 s | 206,829 |
| | 3 | 80.4M | 23.3 GB | 464 s | 133,997 |
| | 4 | 878M (count) | | | 99,241 |
| gnomAD PASS, all (7,584,432) | 2 | 81M (count) | about 25 GB | | 8,517,701 |
| | 3 | 594M (count) | | | 8,039,325 |

Genome-wide GIAB HG001 v4.2.1 (3,893,336 records), measured on October 13, 2026
after the [speed changes](#speed-changes) with
`scripts/benchmark_prefix_hash_scan_vcf.jl ... all` (61 Cas9 guides, distance 3,
eight threads, five alternating runs):

| k | Fragments | Fragment bases | Descriptors | Peak RSS | Time (median, range) | Reference only | Truncated anchors |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 4.42M | 235M | 1.35 GB | 6.5 GB | 49 s (44-100 s) | 2.5 s | 10,052 |

Preparation took a median 45.8 s of the 49 s, and fragment materialization 2.6 s
summed over workers. Output was byte-identical to the genome-wide k = 3 run made
before the speed changes (940 variant rows). The slow runs (84 s and 100 s)
directly followed another VCF run, and their preparation alone took 77-92 s.

Memory is about 0.3 GB per million fragments. Preparation holds the fragments of
every selected chromosome, so a genome-wide run costs roughly the sum over
chromosomes (chr22 is about 1.6% of GRCh38). Run cohort VCFs one chromosome at a
time with `region`, filter by frequency upstream, and use k = 2 for dense cohorts.
For single-sample VCFs, raising k costs little.

Before the cap, exhaustive enumeration of gnomAD AF ≥ 1% chr22 exceeded 177 GB and
was stopped.

## Realistic validation

`test/local_human/vcf_realistic.jl` checks real inputs on GRCh38 with 61 Cas9
guides at distance 3. It needs local data in the ignored
`test/local_human/data/vcf/` directory:

```bash
julia --project=. -t 8 test/local_human/vcf_realistic.jl search  VCF CHROM|all OUTDIR
julia --project=. -t 4 test/local_human/vcf_realistic.jl windows VCF CHROM N SEED [MAXRECORDS]
# VCF_K=k overrides vcf_max_alleles (default: the package default)
```
 Inputs were normalized with bcftools 1.19
(`norm -f GRCh38 -c w -m+any`; gnomAD and GIAB also renamed `chrN` to `N`; ClinVar
records on NT_/NW_ contigs absent from the Ensembl primary assembly were dropped).
No REF mismatches were reported.

- `windows` samples real records in 100-320 bp windows (at most 8-10 records,
  at most 256-1,024 haplotypes), builds guides from ALT-haplotype sites with up to
  d random substitutions, and compares all rows with the independent oracle.
- `search` runs a chromosome or the whole genome and checks that reference rows
  equal a run without VCF, that output is byte-identical across 1/8 threads,
  64 KiB chunks and the portable SIMD backend, and that every variant row matches
  a haplotype rebuilt from its own annotations (sequence, PAM, distance, endpoints,
  touched alleles). Reference-row `overlapping_variants` are recomputed from the VCF.

Results on October 11-12, 2026, all passing:

| Input | k | Windows | Search scope | Variant rows validated |
|---|---:|---:|---|---:|
| GIAB HG001 v4.2.1 | exact | 100 | whole genome | 942 |
| GIAB HG001 v4.2.1 | 3 | 100 | whole genome | 940 |
| ClinVar | exact | 100 | — | — |
| ClinVar | 3 | 100 | chr22 | 1,020 |
| gnomAD AF ≥ 1% | exact | 100 | — | — |
| gnomAD AF ≥ 1% | 3 | 100 | chr22 | 88,071 |
| gnomAD AF ≥ 0.1% | exact / 3 | 100 / 100 | — | — |
| gnomAD all PASS | exact / 3 | 5 / 8 | — | — |

"Exact" windows use a cap equal to the window's record count. All-PASS gnomAD has
about one record per 5 bp, so few windows are small enough for the oracle. At k=3,
GIAB truncated 10,052 anchors genome-wide and lost 2 of 942 variant rows.

## TODO

Backlog toward publication, in suggested order.

Polish:

- [x] Keep `test/local_human/vcf_realistic.jl` in git (data stays ignored).
  `test/variants/regressions.jl` keeps only "reference cap retirement and variant
  deduplication" and "example-genome fixture and CLI". Its other assertions not
  covered elsewhere (region without VCF, legacy backend, counts, preparation
  stats, long-insertion filter/summary agreement) moved to the input, execution
  and corpus tests.
- [x] Throttled progress records for VCF read, enumeration and the variant scan.
- [x] Actionable input errors: `bcftools annotate --rename-chrs` for
  `chr22`/`22` mismatches, `bcftools view -t ^contig` for contigs absent from the
  reference (ClinVar NT_/NW_ records), `bcftools norm -f` for REF mismatches.
- [x] Standalone app (built 2026-10-12): search with `--vcf`, `--region`,
  `--vcf_max_alleles` and `--vcf_annotate`, then `filter` and
  `summarize --variant_aware`, match `CHOPOFF.main` byte for byte with
  `--simd_backend avx2`.
- [x] Standalone app bug, not VCF-specific: with the AVX-512 backend (the `auto`
  choice on a Xeon Gold 6126) the reference scan returned no rows on
  `test/sample_data/vcf/references/contracts.fa`, `--region chr1`, d1. Cause: the
  sysimage code could not inline the AVX-512 profile primitive and passed its
  `<64 x i8>` argument in mismatched registers. The primitive now takes a pointer
  and returns scalars (see prefixHashScan.md, SIMD backend). Rebuilt app
  (2026-10-13): all four backends find the row; GIAB chr22 search, filter and
  `summarize --variant_aware` with `auto` (AVX-512) match `CHOPOFF.main` byte for byte.
- [x] V14: a failure injected into the second batch's variant scan leaves the
  previous output and directory unchanged.
- [x] `align` agrees with a naive edit-distance DP on 5,000 seeded pairs.

Cheap speed and documentation:

- [x] Speed items 1-3 from [Deferred design](#deferred-design), with
  before/after chr22 timings and exact output parity; see [Speed changes](#speed-changes).
- [x] Genome-wide GIAB timing (the main single-sample use case) in the cost table:
  49 s median versus 2.5 s for reference only, after replacing `Base.summarysize`.
- [x] README section: bcftools recipe, one example, output columns, and when to
  change `vcf_max_alleles`. The example runs on the standalone app.
- [x] Line scan for records on unselected chromosomes: GIAB chr22 6.90 s to 2.68 s;
  see [Speed changes](#speed-changes). Standalone app rebuilt 2026-10-13 before
  this change, so it does not include it.
- [x] Profile genome-wide preparation (GIAB, ClinVar): phase times, CPU profile,
  and the slowdown of repeated runs in one process; see
  [Genome-wide preparation profile](#genome-wide-preparation-profile).
- [ ] Read-path fixes 1-2 from the profile (concrete-ALT fast path, allele key
  built once).

Features and evaluation:

- [ ] Genotype and phase support (D01). Restricting walks to alleles carried on
  one phased haplotype gives per-sample results and bounds combinations without an
  arbitrary cap, which makes cohort (1000 Genomes style) runs feasible. Needs its
  own design plan.
- [ ] Head-to-head comparison with CRISPRme and CRISPRitz on 1000 Genomes or GIAB:
  runtime, memory, and found or missed off-targets.
- [ ] Variant-aware count output; VCF searches currently require detail output.
- [ ] Interval/BED regions (D03) and per-chromosome preparation (D07).

## Deferred design

Speed and memory, in priority order. Each needs exact output parity and a measured
win on the chr22 benchmark. Work sizing, stored reference spans, block reads and
the line scan for unselected chromosomes are done ([Speed changes](#speed-changes)).

1. Accept a hit only when its window plus extension overlaps the last taken
   variant of its fragment, and trim fragments to that band. Each site context is
   then verified once instead of once per fragment that contains it.
2. Prepare and scan one chromosome at a time (D07), bounding peak memory by the
   largest chromosome, and share walk prefixes (a trie) instead of copying segment
   paths. These are the changes that could make k = 4-5 feasible for cohorts.

Genotype support needs original ALT indexes, per-record ploidy, missing calls, and
phase-block-aware carrier checks. Intersect haplotype slots within each phase block,
then combine sample-level feasibility across independent blocks. A global slot AND
would reject valid cross-block combinations. Skipped variants hidden by a taken
deletion must not impose reference-allele requirements.

SV support needs junction provenance, mate handling, and combinations with nearby
variants. Generic `<DUP>` does not specify a tandem junction. BND identity must include
orientation and inserted sequence. Short inversions or duplications can place both
junctions inside one site. Reference-only junction flanks would miss combined sites.

Interval/BED selection needs deletion-aware reachability. A fixed W-base reference
margin can omit a distant SNP brought into the interval by a deletion. Whole-chromosome
selection avoids this boundary problem.

The IUPAC compression valve remains deferred. It must account for existing reference
ambiguity and preserve allele provenance. Explicit ALT expansion is the first-release
correctness path. Genotype mode, SVs, sample summaries, `lost`, and retirement of the
legacy vcfDB implementation remain separate work.

## Legacy vcfDB audit (historical)

Source: `src/db_vcf.jl`. Listed by impact.

### False negatives (safety-critical)

1. **REF allele is never enumerated for nearby variants** (`find_ots_many`, l.120-133). Combinations are built from ALT alleles only. Two biallelic SNPs within one site give only ALT1+ALT2; ALT1+REF2 and REF1+ALT2 are lost. `test/sample_data/artificial.vcf` hides this because it lists REF inside ALT (`A -> A,G,T`). Real VCFs do not.
2. **Clustering window too short** (l.288). Grouping uses `l = length_noPAM + distance`, but a site with its PAM spans `lp = length + distance`. Variants `l+1 .. lp-1` bp apart are never combined, which covers PAM-plus-protospacer combinations. Bulges and insertions need further margin.
3. **Multi-ALT IUPAC compression can drop sites** (`parse_vcf`, l.24-32). ALTs of equal length are merged position-wise into IUPAC codes. Every merged base counts toward `motif.ambig_max`, so `findguides` drops the site once a window collects more than `ambig_max` such bases. The merge is also over-inclusive: `AC,GT` becomes `RY`, which also matches `AT` and `GC`.
4. **Long insertions are missed** (`find_ots_one`, l.81). The overlap test checks only the first and last ALT base, so a site lying fully inside an insertion longer than the site is not reported.

### Wrong sequence or coordinates

5. **Overlapping variants are concatenated as if independent** (l.163-171). A deletion covering a SNP produces a negative gap range and a wrong sequence. Overlapping alleles are mutually exclusive.
6. **`ALT = "."` becomes a fake deletion** (l.18-19). `dna""` replaces the REF base. Symbolic alleles (`<DEL>`, `*`, `<NON_REF>`) throw in `LongDNA{4}(...)`.
7. **Reverse-strand hits are reported as `+`** in `find_ots_one` (l.98, `trues(length(g_rev))`).
8. **Positions do not follow the prefixHashDB convention.** prefixHashDB stores `range.stop` for forward Cas9 and `range.start` for reverse (`normalize_to_PAMseqEXT`, `db_helpers.jl:166`).
   - forward in vcfDB: `seq_start + x.stop`, one too high (should be `seq_start + x.stop - 1`)
   - reverse in `_one`: `seq_start + x.start - 1`, correct
   - reverse in `_many`: `seq_start + x.stop`, off by the site length
   - the convention is hard-coded for `extends5 = true`, so Cas12a positions land on the wrong end
9. **No liftover after indels.** Positions downstream of an indel stay in alt-sequence coordinates.
10. **A variant within `lp` of a chromosome edge throws a BoundsError** (`seq_start < 1`).

### Scalability

11. **Clusters are built as chains.** Grouping links variant i to i+1 transitively, so a dense region becomes one cluster and enumeration runs over the whole cluster instead of one site window. With REF included (bug 1 fixed), the cost is 2^n over the chain.

### Search side (shared with prefixHashDB)

12. **Ambiguous-base offset (suspected, needs a test).** `search_vcfDB` l.494 and `db_prefix_hash.jl:445` compute the offset as `sum(is_ambig[1:bv_start]) + 1`; the correct form is `sum(is_ambig[1:bv_start-1]) + 1`. Whenever the first base of a stored site is ambiguous, the decoded bases shift by one. The sum is also O(N) per candidate. This affects plain prefixHashDB too when the reference has IUPAC bases.

### Semantics

13. GT, phasing and FILTER are ignored; every record is assumed present in one hypothetical genome. Output holds only variant sites, so the user must merge them with reference hits.
