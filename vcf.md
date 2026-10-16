# Variant-aware prefixHashScan

`search_prefixHashScan` accepts a normalized VCF through `vcf` or `--vcf`.
It searches sequence-resolved SNVs, MNPs, insertions, deletions,
unequal-length replacements, and precise symbolic `<DEL>`, `<DUP>` and `<INV>`
structural variants (see [Structural variants](#structural-variants)). It searches compatible REF/ALT combinations up to
an ALT-allele cap. With `vcf_mode=:genotype` it searches only combinations that a
sample copy may carry, using GT and PS, and reports carriers. It does not interpret
FILTER or allele frequency.

Status: complete for single-sample VCFs, validated on GIAB HG001 genome-wide.
Cohort and catalog VCFs (gnomAD, ClinVar) work with a small cap. Fragments are
enumerated one chromosome at a time, so union-mode memory follows the largest
chromosome; genotype mode still holds every selected chromosome's GT entries.

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
| V02 | Plain VCF, gzip and BGZF | [x] implemented, reader | tested; [byte parity, genuine BGZF split within a row](test/variants/input.jl), [multithreaded BGZF batches, corrupt and truncated blocks](test/variants/bgzf.jl) |
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
| D01 | Genotypes, phase blocks and carrier/sample summaries | [x] implemented, [genotypes](src/variants/genotypes.jl) | tested; [token property tests](test/variants/genotypes.jl), [per-sample oracle, hand cases, summaries](test/variants/genotype_corpus.jl); see [Genotype mode](#genotype-mode) |
| D02a | Symbolic DEL, DUP and INV and nearby small-variant combinations | [x] implemented, reader/enumerator | tested; [symbolic-literal parity, oracle, skips, errors, genotype mode](test/variants/sv.jl), [90 random junction cases](test/extended/variants.jl), genome-wide HG002 parity; see [Structural variants](#structural-variants) |
| D02b | BNDs, long-range junction anchors, `<INS>` without sequence | [ ] deferred | skipped with a warning |
| D03 | Interval/BED selection with deletion-aware reachability | [x] implemented, [region parser](src/variants/vcf_read.jl), [anchor selection](src/variants/enumerate.jl), scanner | tested; [parser and restrictions](test/variants/input.jl), [filtered-search parity, deletion bridge, early stopping](test/variants/regions.jl); see [Regions](#regions) |
| D04 | IUPAC compression valve | [ ] deferred | untested; explicit expansion only |
| D05 | `lost` summaries | [~] partial | genotype mode reports copies that lose a reference site; union mode unchanged |
| D06 | Legacy vcfDB retirement | [ ] deferred | legacy implementation/tests remain |
| D07 | Parallel, per-chromosome preparation | [~] partial, [enumerator](src/variants/enumerate.jl), [scanner](src/variants/search.jl) | tested; [block parity across threads](test/variants/regions.jl), byte parity on GIAB and ClinVar; the read and GT store still cover all selected chromosomes |
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

SV call sets need the same contig names. Symbolic records need no `norm`; to
search SVs together with small variants, merge them into one sorted VCF
(`bcftools concat -a` then `bcftools sort`).

Keep multiallelic records together when possible. The reader splits them internally
and remembers that their ALT choices are exclusive. Pre-split records without
source provenance become independent choices, subject to footprint conflicts.
That can add combinations absent from the original multiallelic record.

The reader uses only the `##fileformat`, `##INFO`, `##FORMAT` and `#CHROM` header
lines and skips other meta lines, including malformed ones. A tabix index is unnecessary. The reader streams plain VCF or gzipped VCF. BGZF
input (bgzip, `bcftools -Oz`) inflates its blocks on all threads while records
are parsed; other gzip input uses the single-threaded decompressor.
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
comma-separated chromosome names or a BED file path; see [Regions](#regions).
Unknown names and interval strings are errors. Chromosome selection also works
without a VCF, including count output and the legacy scanner.

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

Symbolic SVs follow [Structural variants](#structural-variants). Literal sequence
alleles have no size cutoff. Users must left-align upstream. The reader trims common
padding but does not validate left alignment or merge equivalent records.

Use `vcf_annotate` or `--vcf_annotate INFO/AF,FORMAT/AF` to retain selected fields.
Unknown header fields fail at startup. INFO/FORMAT fields with Number=A select the
chosen ALT value. Number=R retains the REF and chosen ALT values. Other cardinalities
retain their original value. FORMAT annotations name their samples without interpreting GT.

## Structural variants

Precise symbolic deletions, tandem duplications and inversions become literal
alleles of any length. They then combine with nearby small variants, take part in
genotype mode, and are filtered and summarized like other alleles. Every other
SV form is skipped, never an error, and each skipped kind is counted.

| ALT | Searched as | Footprint and sequence |
|---|---|---|
| `<DEL>`, `<DEL:*>` (including `<DEL:ME:*>`) | deletion of POS+1..END | `[POS+1, END+1)`, empty |
| `<DUP:TANDEM>` | insertion of a reference copy before the original | gap POS+1, `ref[POS+1..END]` |
| `<DUP>` | as `<DUP:TANDEM>`, with a warning | gap POS+1, `ref[POS+1..END]` |
| `<INV>` | equal-length replacement | `[POS+1, END+1)`, reverse complement of `ref[POS+1..END]` |
| any of the above with INFO/IMPRECISE | skipped, counted as `imprecise` | breakpoints, and so junction sequence, are uncertain |
| `<INS>`, `<INS:*>` (including mobile elements) | skipped, counted as `ins` | inserted sequence unknown; literal insertions are searched |
| `<CNV>`, `<CNV:*>` | skipped, counted as `cnv` | copy number without arrangement or sequence |
| `<DUP:*>` other than `TANDEM` (e.g. `<DUP:INT>`) | skipped, counted as `dup_other` | location of the copy unknown |
| breakends (`t[p[`, `]p]t`, `.A`, `A.`) and translocations | skipped, counted as `bnd` | needs junctions between distant or reverse-oriented segments (D02b) |
| any other symbol (`<CPX>`, `<CTX>`, ...) | skipped, counted as `other` | complex or undefined arrangement |
| a non-IUPAC letter (`Z`) | error | |

END comes from INFO/END, else POS + |SVLEN| (the ALT's value when SVLEN has one
per ALT). A record without either, or with END not in (POS, chromosome length],
fails. REF is the padding base and is checked like any REF; GATK-SV writes `N`,
which matches every base. `CIPOS`/`CIEND` without `IMPRECISE` are searched at the
stated breakpoints.

IMPRECISE records are skipped on purpose. The target use is well-annotated
genomes of decent quality, and an imprecise breakpoint makes junction sites
guesses. Skipping them also drops what does not depend on the exact breakpoint:
reference sites a carrier loses inside an imprecise deletion are not listed in
`overlapping_variants` or `lost`. The 1000 Genomes SV set has 7,883 such deletions
(plus 1,546 duplications and 471 inversions). Records without the flag are
treated as exact even when they carry `CIPOS`/`CIEND` (about 5,600 in that set). Reference `N` inside a copied or inverted segment stays `N`
and falls under `ambig_max`.

Skips are reported once per run: a warning with the count and up to five example
loci per kind, the counts (`sv_skipped_<kind>`, `sv_del`, `sv_dup`, `sv_inv`) in the
`prefixHashScan variants` log line, and `PrefixHashScanStats.vcf_sv_searched` and
`vcf_sv_skipped`.

The allele key keeps the original ALT and adds END: `chrom:POS:REF><DUP:TANDEM>/END`.
Inserted bases of a duplication have the origin `chrom:POS:REF><DUP:TANDEM>/END:offset`.

Meaning of the rows:

- Every interior site of a carried duplication exists twice. The copy's site is a
  separate variant row on the insertion-origin axis, so `summarize --variant_aware`
  counts both copies. Offsets within `distance` of a copy end join the reference
  axis at POS, as for literal insertions.
- The new site at the duplication junction (`ref[END]` followed by `ref[POS+1]`) and
  sites bridging a deletion are variant rows like those of literal alleles.
- An inversion is reported per base, as for MNPs. Reference rows inside it list it
  in `overlapping_variants`, and in genotype mode carriers lose them; the inverted
  site appears as a variant row.

Limitations:

- Small variants inside a duplication apply to the original copy only; the inserted
  copy is reference sequence.
- Small variants inside an inversion conflict with it, so their combinations inside
  the inverted segment are not searched.
- Reference sites crossing gap POS|POS+1 of a duplication list it in
  `overlapping_variants` and `lost`, although the tandem copy keeps that adjacency.
- Union-mode summaries count overlapping duplication alleles of one region
  (different END) separately, as for alternative insertions. Use
  `summarize --per_sample` for cohort SV VCFs.
- Interior rows of large duplications spend the variant early-stopping budget.
- Breakends and translocations need junctions between distant or reverse-oriented
  segments (D02b) and are skipped.

### SV validation

[test/variants/sv.jl](test/variants/sv.jl) runs each symbolic case twice: as the
literal VCF, whose rows the [oracle](test/helpers/vcf.jl) checks, and as the symbolic
VCF. The rows must be equal once allele keys are renamed. Cases: a 12 kb deletion
bridging a site with and without a flank SNP, tandem duplications and an inversion
longer than the scan chunk, an upstream walk taking a long duplication, deep copy
sites in filter and summary, SVLEN-only and multiallelic records, every skipped kind,
malformed records, and genotype mode. The extended suite adds 90 seeded random
cases (seed `20261009`) in which a DEL, DUP or INV joins the two halves of a site at
a random split, with up to three SNPs and a guide with up to `d` substitutions;
every case reports a row using the SV.

Real data, GRCh38, 61 Cas9 guides, distance 3, eight threads, measured
October 15, 2026 on a loaded host (files and preparation in `../sv_data/README.md`):

- GIAB HG002 v5.0q, small variants and SVs, all literal (5.9M records after
  normalization). A copy with 12,430 deletions of at least 50 bp written as `<DEL>`
  and 788 insertions that repeat the following reference written as `<DUP:TANDEM>`
  gave the same 28,158 rows genome-wide as the literal file, once keys were renamed;
  23 rows use a converted SV. Both runs built 124.4M fragments (6.87 Gb), truncated
  452,116 anchors at k = 3, and took 269 s (literal, first) and 214 s, with 14-16 GB
  peak RSS. This v5.0q file is far denser than the v4.2.1 HG001 set above.
- 1000 Genomes high-coverage SV set, chromosome 22 (2,643 records: 1,449 `<DEL>`,
  14 `<DEL:ME*>`, 512 `<DUP>` up to 3.4 Mb, 4 `<INV>`; IMPRECISE records and the
  other kinds skipped). Union mode, `vcf_realistic.jl search`: 183 variant rows, all
  second copies of reference sites inside duplications, each rebuilt from its
  alleles; reference rows equal the search without VCF; `overlapping_variants`
  correct on all 457 reference rows; output byte-identical across 1/8 threads,
  64 KiB chunks and the portable backend; 8.9 s. Genotype mode with all 3,202
  samples, `gtsearch`: 18,786 certain carrier copies and 11,001 lost reference
  copies, all consistent with GT; 13.7 s, 1.8 GB peak RSS.
- The 1000 Genomes header contains `###INFO=<ID=STRANDS,...>`, which bcftools
  accepts and the Julia VCF header parser rejects. The reader now passes only
  `##fileformat`, `##INFO`, `##FORMAT` and `#CHROM` lines to that parser and skips
  every other meta line.

Existing outputs did not change: GIAB HG001 and ClinVar chromosome 22 searches
were byte-identical before and after the walk changes.

Long alleles stay cheap: a walk passes a taken deletion in one step, and a taken
allele that runs past the walk's goal is cut there, because no site overlapping
the anchor reaches beyond it. Each duplication or inversion is scanned in full
once, as its own anchor. `ScanVariant` stores each copied or inverted sequence, so
memory grows with their total length.

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
to a concrete ALT. Empty ALT segments retain deletion junctions. Preparation reads
the VCF once and keeps the variants for all guide batches. The variant scan then
enumerates fragments one chromosome at a time, in blocks of consecutive anchors
(four per scan thread) on all scan threads, scans them and frees them. Each block
reads the merged reference spans its fragments use once, through a forward 1 MiB
block reader. Workers materialize fragments from these spans into reusable
buffers without file reads, then use the existing prefix kernel, Myers verifier,
and traceback. Enumeration repeats for every 64-guide batch, so a search with more
than 64 guides pays it once per batch.

Variant work packs each chromosome's fragments into items of about a quarter of its
fragment bases per scan thread (at least 64 KiB, at most the chunk size, 2 MiB by default), with
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

## Genotype mode

`vcf_mode=:genotype` (CLI `--vcf_mode genotype`) reads GT and PS for all samples or
for `vcf_samples` (`--vcf_samples S1,S2`). A walk takes an allele only when some
copy of a selected sample may carry it together with every allele already taken,
so the walk follows carried haplotypes. The cap defaults to unlimited in this mode.
ALTs that no selected sample may carry are dropped, so they appear in neither
walks nor `overlapping_variants`.

Phase groups:

- `|` without PS: one phase set per sample and chromosome (reference-panel phasing).
- `|` with PS: one phase set per PS value.
- `/`, partial phasing (`0|1/2`) and haploid calls: their own group, copies permute freely.
- VCF 4.4 leading indicators (`|0|1`) are read; PSL and PSQ are ignored.

A copy carries a variant row when it carries every allele the site uses and none
of the other alleles that intersect a consumed reference base or insert at a gap
the site crosses (including the edges of a taken deletion). Alleles hidden inside
a taken deletion impose nothing. A copy loses a reference row when it carries an
allele in `overlapping_variants`.

Genotype output adds six columns:

| Column | Meaning |
|---|---|
| `carriers` | Variant rows: copies carrying the site |
| `lost` | Reference rows: copies lacking the site |
| `n_carriers`, `n_lost` | Certain copies (`S:j` and `S:*` tokens) |
| `n_possible`, `n_lost_possible` | Samples with an `S:?` token |

Tokens consider every permutation of phase groups and every value of missing
alleles (`.`):

- `S:j` holds for copy j in every case. `S:j@PS` names copy j of a PS block;
  without a suffix j is the chromosome-wide copy. Indexes of different PS blocks
  are not comparable.
- `S:*` holds for some copy in every case, but which copy is unknown, for example
  a single unphased heterozygous allele. It repeats once per certain copy.
- `S:?` holds in some case only: unphased combinations, alleles in different PS
  blocks, missing alleles. In `lost`, `S:?` after certain tokens means further
  copies may lose the site.

Rows carried by no selected copy are dropped. Some call sets, including the 1000
Genomes panel, place overlapping alleles on one copy (for example two deletions of
one STR) instead of using `*`. The run warns once with the pair count. Walks never
take both alleles, and by the overlap rule only rows of the allele that hides the
other credit that copy. For ploidy above two, `S:?` can
appear where a counting argument proves some copy carries the site.

`summarize_offtargets(res; per_sample=true, vcf, vcf_samples, distance)` (CLI
`summarize --per_sample --vcf in.vcf.gz`) counts loci per guide and sample at the
best distance among rows the sample carries. `D` columns count loci some copy
certainly carries. `P` columns are the upper bound with possible carriers. A
reference row is possible while fewer than ploidy copies certainly lost it, and
certain when `lost` also has no `S:?` for the sample. The VCF supplies sample names
and per-chromosome ploidy, since hom-ref samples appear in no token.

### Genotype validation and cost

[test/variants/genotypes.jl](test/variants/genotypes.jl) checks the GT parser,
the store and 3,000 random token cases against brute-force enumeration of phase
permutations and missing alleles. [test/variants/genotype_corpus.jl](test/variants/genotype_corpus.jl)
runs hand cases (cis, trans, unphased, PS blocks, deletion edges, sample
selection, conflicts, errors, CLI) and random cases against a
[per-copy oracle](test/helpers/vcf_genotypes.jl). The oracle builds each copy's
haplotype per assignment, searches it with linearDB (see
[Validation and measurement](#validation-and-measurement)) and derives tokens from
the definitions. Rows must also be a subset of union-mode rows.
Random cases default to 20 (`CHOPOFF_VCF_GT_CASES`, `CHOPOFF_VCF_GT_SEED`); the
extended suite runs 200. On 2026-10-14, 200 cases with seed `20261014` passed.

Real data checks use `test/local_human/vcf_realistic.jl gtsearch|gtwindows`, which
rebuild each listed copy's haplotype from the VCF and check the reported site:

- NA12878 from the 1000 Genomes 30x phased panel, chr22 (59,328 records): 12
  variant rows with 16 certain carriers, 457 reference rows with 25 lost copies,
  no failures. A copy with synthetic 50 kb PS blocks and 10% unphased hets also
  passed. GIAB HG001/HG002 v4.2.1 are unphased (`/`, no PS) and were not used.
- 1000 Genomes chr22, every 32nd of 3,202 samples (100): 127 variant rows, 2,698
  certain carriers, 2,895 lost copies, no failures. 100 random windows with 20
  samples also passed.

Reference rows matched the reference-only search, and outputs were byte-identical
across one or eight threads and 64 kb chunks.

Cost on GRCh38 chr22 (61 Cas9 guides, distance 3, eight threads, no allele cap,
shared host under load). Reference alone took 7.4 s.

| Samples | Wall | Peak RSS | Fragments | Max per anchor | GT entries | Descriptors |
|---:|---:|---:|---:|---:|---:|---:|
| 1 (NA12878 VCF) | 8.5 s | 1.2 GB | | | | |
| 100 | 72, 91 s | 1.5 GB | 0.55 M | 6,772 | 6.3 M | 0.28 GB |
| 3,202 | 225, 253 s | 6.0 GB | 1.75 M | 8,192 | 205 M | 4.2 GB |

Each run is a fresh process including compilation. Reading inflates BGZF blocks
and parses GT columns on all threads, in batches of 1,024 lines; the rest of each
record stays sequential. Alternating runs against sequential GT parsing gave 111
and 115 s (100 samples) and 357 and 422 s (3,202 samples); before that, the
single-threaded gzip reader added another 5-15%. Output was byte-identical in all
runs. Reading the whole file line by line takes 7.7 s with the BGZF reader and
29 s with the single-threaded one. At 3,202 samples, garbage collection takes
50-65 s, and the main thread spends about 70 s reading lines and checking records.
The variant scan took under 10 s.
Union mode needs the cap on the same file; see
[Allele cap and measured cost](#allele-cap-and-measured-cost).

Genotype mode keeps the GT store for every selected chromosome for the whole
run. The store is sparse (205M entries for 3,202 samples on chr22, not samples
times records), but scaled by genome share a genome-wide 1000 Genomes run would
need about 13 billion entries and roughly 250 GB. Run cohort genotype searches
one chromosome at a time with `region`; the 1000 Genomes panel ships per
chromosome anyway. Single-sample genome-wide runs hold about one entry per record.

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

## Regions

`region` (CLI `--region`) also accepts a BED file path. The reader skips `#`,
`track` and `browser` lines, reads the first three tab-separated columns
(0-based, half-open), and merges overlapping or adjacent intervals. Unknown
chromosomes, `start >= end` and intervals past the chromosome end are errors.
BED input must be uncompressed.

A row is reported when its site overlaps an interval. The site runs from `start`
to `end`, including the PAM. Variant rows use their lifted endpoints, so a site
that bridges a deletion covers the deleted reference between its endpoints, and
an endpoint inside an insertion counts at POS. Reference-row
`overlapping_variants` and genotype `lost` still consider every variant on the
chromosome.

The reference scan reads only candidate windows that can reach an interval: each
interval widened by the candidate window and `distance`. Verification checks the
exact traceback span before a row spends early-stopping budget, so rows outside
the region never retire a guide.

Variant anchors are selected with a deletion-aware bound. After its anchor a walk
emits at most W bases. Reference bases outside every variant footprint can only
be passed by emitting them, so the walk ends before the (W + 1)-th such base.
Deletions therefore extend the reach without a fixed margin, and the bound
assumes no allele cap, so it also holds in genotype mode. Records on the selected
chromosomes are read whole; only the anchors are restricted.

BED regions need detail output and a streaming backend. Count output has no
traceback span to test, and the legacy scanner has no interval support; both are
errors.

## Validation and measurement

The regular entrypoint is [test/src/variants.jl](test/src/variants.jl). Its
[oracle](test/helpers/vcf.jl) constructs complete alternative haplotypes with a
separate edit model, searches them with linearDB, and lifts coordinates with its
own per-base map. linearDB is the original CHOPOFF engine and shares no search code
with prefixHashScan: site finding, filtering, alignment and offsets are all separate.
It compares every output column except the two alignment strings, including
distance, allele annotations, reference overlaps and insertion origins, and rejects
duplicate rows. Alignment strings are skipped because linearDB (`suffix_align`) and
prefixHashScan (`align`) may place a gap differently in an alignment with the same
distance. Haplotype sets with any IUPAC base use the prefixHashScan brute-force
backend instead: linearDB counts ambiguous bases over guide and PAM only and lets an
`N` in the extension match freely, while prefixHashScan counts the whole window. A
parity run on 2026-10-16 over the whole variant suite found the same sites and
distances in both engines.
Hand-calculated checks validate the oracle's edit model. Production VCF parsing,
enumeration, conflict checks and coordinate lifting do not participate in the oracle.
[test/src/distance_metrics.jl](test/src/distance_metrics.jl) checks `align` against
a naive edit-distance DP on 5,000 seeded guide/reference pairs with k = 1-4,
including IUPAC reference bases and indels; distances, gapless sequences and
column costs all agree.

[Permanent fixtures](test/sample_data/vcf/README.md) combine the existing
semirandom genome example with tiny companion references. Every positive sequence
case requires a specific allele witness, so an empty result cannot pass. Invalid
inputs are separate files. The generator's `--check` mode regenerates into a
temporary directory and checks bytes and file lists without changing the corpus.

[Real-data fixtures](test/sample_data/vcf/real/README.md) are 98 GRCh38 windows
(about 19 kb, under 1 MB with expected outputs) cut from GIAB HG002 v5.0q, GIAB
HG001, the 1kGP phased panel and PS file, the 1kGP SV set (trio
NA12878/NA12891/NA12892), ClinVar and gnomAD. `make_real_vcfs.jl` picks windows
greedily until each tag (variant kind, size class, GT pattern, PS switch, SV kind,
overlapping, adjacent and clustered records) is covered twice, and needs the large
local files. [test/variants/real.jl](test/variants/real.jl) compares every window
with the haplotype oracle, checks genotype carriers with the per-copy oracle where
the window allows it, checks HG002 symbolic SVs against their literal twin, and
compares whole-set union, genotype and `summarize --per_sample` output with
`expected/`. It runs in about 100 s.

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

After the region and per-chromosome changes (2026-10-15, four threads), the
regular suite passed 10,942 assertions, the extended suite 1,195 plus 231 genotype
assertions (200 cases), fixture regeneration matched, and the full package suite
passed.

The annotation regression exposed cumulative ID ranges in the installed VCF parser.
The reader preserves the complete ID field so multiple IDs appear exactly once in
both variant and reference-overlap annotations.

`PrefixHashScanStats` records input counts, fragment counts and bases, the largest
fragment count of one anchor (`max_anchor_fragments`), anchors truncated by the cap
(`vcf_truncated_anchors`), descriptor bytes (variants plus the largest
chromosome's fragments and reference spans), preparation time (read plus
enumeration), variant fragment materialization time (`variant_io_ns`), candidates,
anchor rejections, accepted variant rows, and retired variant guides.
`verbose=true` reports read counts once before guide batching and fragment counts
after the first batch's variant scan. It also logs throttled
`prefixHashScan progress` records (every `progress_interval` seconds) with
`batch` set to `VCF read`, `<batch> VCF enumerate` or `<batch> VCF` for the variant scan.
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

Genome-wide ClinVar at k = 3 was not practical: 455M fragments, 91 GB of
descriptors and 158 GB peak RSS, with enumeration taking 95% of the time. After
[per-chromosome enumeration](#per-chromosome-enumeration) the full search takes
460 s and 29 GB.

Repeated GIAB preparations in one process took 40.6-44.0 s, so the 84 s and 100 s
benchmark runs were not reproduced with preparation alone. Their cause (host load,
or GC state left by the previous search) is unknown.

Applied after the profile, measured October 13, 2026 (same benchmark, five
alternating runs on chromosome 22, three genome-wide):

- Plain ACGT ALTs skip the IUPAC expansion in `foreach_concrete_variant`.
- The `chrom:pos:ref>alt` key is built once per ALT; `variant_token` skips
  `replace` when the string has no `%`, `;` or `|`.
- ID, REF and ALT come from the tab-split line. `VCF.Record` runs only when
  `vcf_annotate` needs INFO or FORMAT fields, so without it selected records no
  longer get the package's QUAL, FILTER and INFO syntax checks.
- CHROM is a `SubString`; a `String` is allocated only on a chromosome change.
- The enumeration walk is a top-level function. As a recursive local closure it
  captured itself and `capped` in `Core.Box`, so every call dispatched
  dynamically.
- Reference span bytes are reserved before reading.

| Input | Before | After | Preparation before / after | Allocated before / after |
|---|---:|---:|---:|---:|
| Synthetic chr22 | 0.99 s | 0.86 s | 0.41 / 0.21 s | |
| GIAB chr22 | 2.86 s | 2.64 s | 2.26 / 2.03 s | |
| ClinVar chr22, k = 3 | 28.0 s | 20.8 s | 26.0 / 17.7 s | 11.1 / 6.8 GB |
| GIAB genome-wide | 55.9 s | 32.0 s | 50.4 / 26.7 s | 19.6 / 10.1 GB |

Output was byte-identical for all four, and prepared data (variants, counts,
fragments, segments, reference spans) was identical for every valid fixture,
including `vcf_annotate`.

### Per-chromosome enumeration

Measured October 15, 2026 with `scripts/benchmark_prefix_hash_scan_vcf.jl` (61 Cas9
guides, distance 3, eight threads, k = 3), one measured run after warmup per
process; genome-wide GIAB alternated three processes per version. The host was
loaded, so chr22 times are higher than in the tables above. Variant output was
byte-identical to the previous code in every run.

| Input | Before | After | Preparation before / after | Peak RSS before / after |
|---|---:|---:|---:|---:|
| GIAB chr22 | 6.8 s | 6.3 s | 6.2 / 5.6 s | 1.34 / 1.31 GB |
| ClinVar chr22 | 22.4 s | 12.9 s | 20.4 / 11.2 s | 5.9 / 5.9 GB |
| GIAB genome-wide | 49.6 s (48.7-52.2) | 44.2 s | 46.0 / 38.7 s | 5.2 / 4.0 GB |
| ClinVar genome-wide | preparation alone 1,034-1,386 s | 460 s | | 158 / 29.4 GB |

Preparation here is the read plus enumeration. Descriptors for genome-wide GIAB
fell from 1.38 GB to 0.67 GB (variants plus the largest chromosome's fragments).
Genome-wide ClinVar holds 14.7 GB of descriptors and reports 208,607 variant rows.
Concatenating the thread blocks raised ClinVar chr22 peak RSS to 8.6 GB, so the
scan reads blocks directly. Without a collection between chromosomes, genome-wide
GIAB peaked at 8.6 GB; an incremental collection removed it without a measurable
time cost, and a full one added 18 s.

BED timing, GIAB genome-wide, 1,000 random 1 kb intervals on chromosomes 1-22
and X, three in-process repetitions: reference-only 0.48-0.64 s, with VCF 39.5-40.0 s
versus 38.3-38.4 s for the whole genome. The BED output equalled the whole-genome
rows filtered by site overlap. With a VCF the read of every record on the selected
chromosomes dominates; restricting the read to records near intervals is not
implemented.

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
After the [preparation fixes](#genome-wide-preparation-profile) the median is
32.0 s (31-44 s over three runs), with 26.7 s of preparation and 6.5 GB peak RSS.

Memory is about 0.3 GB per million fragments. Fragments exist one chromosome at a
time, so a genome-wide run costs the largest chromosome's fragments plus the
variants of every selected chromosome; time still sums over chromosomes (chr22 is
about 1.6% of GRCh38). Filter cohort VCFs by frequency upstream and use k = 2 for
dense cohorts. For single-sample VCFs, raising k costs little.

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
  see [Speed changes](#speed-changes). The standalone app rebuilt on 2026-10-13
  (21:44) includes this change and the preparation fixes below; the full package
  suite passed in its build.
- [x] Profile genome-wide preparation (GIAB, ClinVar): phase times, CPU profile,
  and the slowdown of repeated runs in one process; see
  [Genome-wide preparation profile](#genome-wide-preparation-profile).
- [x] Read-path and enumeration fixes from the profile: GIAB genome-wide 55.9 s
  to 32.0 s, ClinVar chr22 28.0 s to 20.8 s; see
  [Genome-wide preparation profile](#genome-wide-preparation-profile).

Features and evaluation:

- [x] Genotype and phase support (D01); see [Genotype mode](#genotype-mode).
- [ ] Head-to-head comparison with CRISPRme and CRISPRitz on 1000 Genomes or GIAB:
  runtime, memory, and found or missed off-targets.
- [ ] Variant-aware count output; VCF searches currently require detail output.
- [x] Interval/BED regions (D03) and per-chromosome enumeration (D07); see
  [Regions](#regions) and [Per-chromosome enumeration](#per-chromosome-enumeration).

## Deferred design

Speed and memory, in priority order. Each needs exact output parity and a measured
win on the chr22 benchmark. Work sizing, stored reference spans, block reads and
the line scan for unselected chromosomes are done ([Speed changes](#speed-changes)).

1. Accept a hit only when its window plus extension overlaps the last taken
   variant of its fragment, and trim fragments to that band. Each site context is
   then verified once instead of once per fragment that contains it.
2. Share walk prefixes (a trie) instead of copying segment paths. With
   per-chromosome enumeration (done), this could make k = 4-5 feasible for cohorts.
3. Stream the VCF one chromosome block at a time, with all guide batches per
   chromosome, so the variants and the GT store are bounded too. This matters for
   genome-wide cohort genotype runs; see [Genotype mode](#genotype-mode).

Genotype support needs original ALT indexes, per-record ploidy, missing calls, and
phase-block-aware carrier checks. Intersect haplotype slots within each phase block,
then combine sample-level feasibility across independent blocks. A global slot AND
would reject valid cross-block combinations. Skipped variants hidden by a taken
deletion must not impose reference-allele requirements.

Precise DEL, DUP and INV are searched (D02a). Breakends remain (D02b): they need
fragment segments that carry their own chromosome and orientation, a two-sided walk
outward from each junction so small variants on both sides combine with it, mate
pairing, and `end_chromosome`/`end_strand` output columns. Long duplications and
inversions could then also be searched by their junctions instead of full copies.

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
