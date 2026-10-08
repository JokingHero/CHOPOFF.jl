# Variant-aware off-target search (VCF)

Status: analysis and design. No code exists yet.

Primary use case: a sequenced real cell profile, human or mouse, supplied as a VCF, searched for every possible off-target site that cell could contain. The design is optimized for that case.

Also supported, without extra machinery: population VCFs (gnomAD, 1000G, HGDP), phased or unphased calls, somatic and mosaic calls, and non-diploid organisms.

Target engine: `search_prefixHashScan` only. `vcfDB` will be retired; section 1 is kept as a list of lessons.

Goal: sensitivity first. Every off-target site that a carrier genome can contain must be found, including combinations of SNVs and indels inside one site. Extra runtime is acceptable, and so are false positives. False negatives are not.

## 1. Current vcfDB: bugs found

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

## 2. How other tools handle variants

**CRISPRme / CRISPRitz** (Pinello lab, Nat Genet 2023):
- SNVs: an "enriched genome" with IUPAC codes (G/A becomes R), searched once.
- Indels: each indel becomes its own small "fake chromosome" with ±25 nt flank, searched separately.
- After the search, IUPAC hits are expanded to concrete alleles and checked against samples through a variant-to-sample hash table. Phased VCFs (1000G) keep only observed haplotypes. Unphased VCFs (HGDP) are checked per sample but genotypes cannot be combined into haplotypes. For gnomAD, populations are treated as pseudo-samples.
- Stated limitation: no structural variants, and no "combinations of indels and SNPs". An indel and a nearby SNP are never enumerated together.

**CRISPR-HAWK** (Pinello / InfOmics, bioRxiv Dec 2025):
- Phased: builds each individual's haplotype sequences and inserts SNVs and indels together.
- Unphased: encodes positions as IUPAC first, then enumerates all feasible haplotypes, but only at candidate guide sites. Equivalent haplotypes are hash-deduplicated across 76k gnomAD genomes.
- Focus is **on-target** guide design. For off-targets it calls CRISPRitz or writes CRISPRme input, so it inherits the indel+SNP limitation.

Gap: neither tool enumerates SNV+indel combinations exhaustively at off-target sites with genotypes optional. The design below covers that case.

## 3. VCF input contract

CHOPOFF does not left-align, merge or filter VCFs. It validates them, repairs only the two trivial cases below (multi-ALT splitting and IUPAC expansion), and fails with a clear message otherwise. Users run the preprocessing themselves, which keeps the parser small and keeps the well-tested tools in charge of the hard parts.

### Required preprocessing

```bash
# 1. split multi-allelic records, left-align and trim indels
bcftools norm -f reference.fa -m -any -c w input.vcf.gz -Oz -o norm.vcf.gz

# 2. sort and index
bcftools sort norm.vcf.gz -Oz -o sorted.vcf.gz
bcftools index -t sorted.vcf.gz

# 3. select which variants are in play (all optional, all the user's choice)
bcftools view -f PASS sorted.vcf.gz -Oz -o final.vcf.gz          # passing records only
bcftools view -i 'INFO/AF>0.001' final.vcf.gz -Oz -o common.vcf.gz
bcftools view -s PATIENT1 cohort.vcf.gz -Oz -o patient1.vcf.gz   # one sample
bcftools view -r chr22 final.vcf.gz -Oz -o chr22.vcf.gz          # one region
bcftools view -R targets.bed final.vcf.gz -Oz -o targets.vcf.gz  # a BED mask
```

Step 3 is where frequency and sample selection happen. CHOPOFF searches whatever survives it and never filters on its own (section 7). `bcftools` is the tool for this; `bedtools` enters only to build the BED mask that `-R` consumes.

Step 1 matters most. `-m -any` splits multi-allelic records into one ALT per line, which is exactly the one-variant-per-record model in section 4, and it removes the IUPAC compression that causes vcfDB bug 3. `-c w` warns on REF mismatches instead of failing.

### Contract that CHOPOFF checks

| Requirement | Check | On violation |
|---|---|---|
| Coordinate-sorted per chromosome | monotonic POS, and each chromosome in one contiguous block | error |
| Chromosome names match the reference | lookup in `dbi.gi.chrom` | skip that chromosome; **error if no chromosome matches** |
| REF matches the reference sequence | compare REF to the genome | skip that record; **error if more than 1% and at least 10 records mismatch** |
| One ALT per record | `count(',', ALT) == 0` | split internally, count it |
| No IUPAC in REF or ALT | alphabet check | expand internally, count it |
| Indels left-aligned | not checkable cheaply | assumed; see below |

The two `error` rows are there because the alternative is worse than failing. A `chr1` versus `1` naming mismatch, or a VCF called against a different reference build, makes every record skip, and the run then finishes "successfully" having searched no variants at all. A clean exit with an empty variant set is indistinguishable from a genome with no variant off-targets, so it must be an error.

Left-alignment is assumed rather than checked. A non-left-aligned indel still yields the correct sequence, but it sits at a different POS, so it may be grouped with a different set of neighbours and may fail to merge with the same indel written the canonical way. The result is a correct site reported at a shifted coordinate, and in the worst case a missed combination with an adjacent variant. `bcftools norm` removes the problem, which is why it heads the recipe.

Every warning is counted and reported once with its total, never printed per record. A VCF with a million unsplit multi-allelic lines must not produce a million warning lines.

Records CHOPOFF accepts but ignores, each counted in the run summary:
- `ALT = "."` or `ALT = "*"` (no sequence; vcfDB turned these into fake deletions, bug 6)
- symbolic `<INS>` and `<INS:ME:*>` without a sequence, and `<CNV>`. `<DEL>`, `<DUP>`, `<DUP:TANDEM>` and `<INV>` are supported, because `END` or `SVLEN` makes their sequence resolvable (section 5)
- any record with CIPOS or CIEND wider than 0 (see section 5)
- `<NON_REF>` from gVCF: a gVCF must be passed through `bcftools view -e 'ALT="<NON_REF>"'` first

### Source formats and what each one needs

| Source | Genotypes | Mode | Notes |
|---|---|---|---|
| 1000G phased | phased GT, `\|` | `:genotype` phased | per-sample haplotypes available; `PS` not needed, phasing is chromosome-wide |
| HGDP, in-house WGS | unphased GT, `/` | `:genotype` unphased | combinations limited to alleles the sample actually calls |
| gnomAD | none (AF only) | `:all` | the only correct mode; no sample set exists, so `:genotype` is a startup error. Name the AF fields in `--vcf_annotate` |
| Somatic (Mutect2, Strelka) | GT present but meaningless | `:all` | VAF lives in `FORMAT/AF` or `AD`; annotate it, never filter on it |
| Long-read phased (WhatsHap, HiPhase) | phased with `PS` blocks | `:genotype` phased | phasing is only valid inside a `PS` block. Two variants in different blocks must be treated as unphased |
| gnomAD converted to pseudo-samples | fake GT, one per population | `:all` | `:genotype` runs but loses real haplotypes (section 7). Treat the file as AF-only |
| Plant or polyploid | GT with more than 2 alleles | either | works unchanged; arity is read per record, no option to set |
| Clinical panel or small region | either | either | the common case; cost is a non-issue |

`PS` block handling is the one genuinely tricky case: a phased VCF is only phased within a block, so `:genotype` phased mode must fall back to the unphased rule for a pair of calls whose `PS` values differ. Phase is decided per sample and per record: a phased file can still contain `/` calls, and those pairs use the unphased rule too (section 4, Genotype modes). Getting this wrong loses real haplotypes.

### Fields read

Nothing else is parsed, which keeps CHOPOFF tolerant of the heavy INFO fields that break vcfDB today.

- always: `CHROM`, `POS`, `ID`, `REF`, `ALT`
- `INFO`: `SVTYPE`, `END`, `SVLEN`, `MATEID`, `CIPOS`, `CIEND` for SVs; plus any key named in `--vcf_annotate`
- `FORMAT`: `GT` in `:genotype` mode; plus any key named in `--vcf_annotate`
- `FILTER` is never read. Use `bcftools view -f PASS` upstream
## 4. Design: local walk enumeration

### Model

The reference is the backbone and each VCF record (one ALT, after splitting) is a **variant**: an on/off choice. A walk is the reference with some set of non-conflicting variants applied.

Each variant has a reference footprint, with the shared padding base stripped:
- a **substitution** (SNP or MNP) or **deletion** occupies the half-open base interval `[s, e)` it replaces;
- an **insertion** occupies the gap before base `s`.

Two variants **conflict** when their base intervals intersect, when one's gap lies strictly inside the other's base interval, or when two insertions share a gap. Two ALTs split from one multi-allelic record always conflict through the same rule. Conflicting variants never appear in the same walk. Everything else combines: a SNP at `p` and an insertion after `p` stay independent, which keeps the real haplotype `C` followed by `GG` from `REF=A ALT=C` and `REF=A ALT=AGG` at `p`.

There are no merged multi-allele groups. Conflict is not transitive: with a deletion of 10-20, a SNP at 12 and a deletion of 18-25, the SNP and the second deletion conflict with the first deletion but not with each other. Merging the three into one exclusive group would drop the SNP plus second deletion walk.

Variants are ordered by walk order: `(s, kind)`, with a base footprint before a gap footprint at the same `s`, then by record order. So a SNP at `p` comes before an insertion after `p`.

The **anchor span** of a variant is the part of the walk that differs from the reference:
- substitution: the substituted bases;
- insertion: the inserted bases, without the padding base;
- deletion: the two bases either side of the join.

Every off-target site that uses at least one ALT allele is a window of the walk of length at most `W = motif span including PAM + distance`, and it overlaps the anchor span of at least one taken variant.

### Enumeration rule: leftmost-ALT anchor

Each such site has exactly one leftmost taken variant whose anchor span it overlaps. Anchoring on it gives each site exactly one owner.

For every variant `i`:
- left flank: `W` reference bases before the anchor span of `i`. Variants there are not taken.
- right side: a depth-first walk over the variants after `i`. At each variant, branch on skip, and on take if it conflicts with nothing already taken. A branch ends when the walk extends `W` bases past the end of the anchor span of `i`. Lengths are measured along the walk, so the set of reachable variants depends on the earlier choices.
- emit one fragment per leaf.

```
dfs(walk, next, taken):
    if walk extends >= W past anchor(i): emit(walk); return
    v = next variant in walk order
    dfs(walk + ref up to v.end, v+1, taken)                      # skip
    conflicts(v, taken) || dfs(walk + alt(v), v+1, taken ∪ {v})  # take
```

The cost per anchor is the number of leaves, `2^m` for `m` mutually compatible variants within `W`. The exponent is the number of variants inside one window, never the size of a dense region.

Two invariants follow from a flank of `W` on each side, and both are worth asserting in tests:
- **Completeness.** A site overlapping the anchor span ends at most `W - 1` bases past it and starts at most `W - 1` bases before it, so every such site fits inside the fragment.
- **No separator contact.** By the same bound, a candidate window that overlaps the anchor cannot reach a fragment edge.

### Anchor test

A hit belongs to its fragment when its **candidate window** (the motif window plus the `distance` extension that Myers reads) intersects the anchor span. This is a range intersection that runs before Myers, so rejected windows cost nothing more and the test works without traceback.

The candidate window is wider than the aligned span. A hit whose best alignment avoids the anchor can therefore also be emitted by the reference scan or by another anchor. That is an extra row, never a lost one, and the commit-time dedup on complete records removes exact duplicates.

### Pseudocode

```
# preparation
variants = read VCF, split multi-ALT, strip padding, compute footprint and anchor span
sort variants by walk order

# enumeration, once per run, reused by every guide batch
for chromosome, for variant i:
    for walk in dfs(i):
        mode == :genotype && !carried(taken(walk), samples) && continue
        push!(descriptors, Fragment(segments(walk), anchor = span of i))
cut descriptors into variant work items of about 2 MiB

# search
genome work items, then variant work items, through the existing scheduler
for variant item:
    build bytes from RefSeg range reads and LitSeg literals, separator between fragments
    run the existing kernel
    for candidate:
        frag = fragment_of(candidate.start)          # searchsortedlast
        intersects(candidate.window, frag.anchor) || continue
        Myers, early-stop counter, traceback
        (chrom, start, end, strand) = lift(frag, aligned span)
        mode == :genotype && (carriers = carriers(span alleles)) is empty && continue
        emit row with the taken variants inside the aligned span
```

### Fragment representation

```julia
struct RefSeg; chrom::UInt; start::Int; stop::Int; isplus::Bool; end
struct LitSeg; seq::LongDNA{4}; variant_id::Int; end

struct Fragment
    segments::Vector{Union{RefSeg, LitSeg}}
    anchor::UnitRange{Int}    # offset of the anchor span inside the fragment
end
```

Fragments are stored as **descriptors**, not sequence, in flat arrays: one `segments` array plus a per-fragment offset range (the layout `is_ambig` uses in `db_prefix_hash.jl`). A patient WGS run makes 5-10M fragments, and one heap-allocated `segments` vector each would dominate memory and allocation time. The struct above is the conceptual view that `fragments.jl` exposes.

Descriptors are built once per run, before guide batching, and reused by every batch. Enumeration parallelizes per chromosome and is merged in chromosome order, so descriptor order, and therefore output order, is deterministic. Enumeration needs no reference sequence: lengths come from the variants. REF validation (section 3) reads reference ranges through the same FAI/2bit handle the scan uses.

`RefSeg` carries a chromosome and a strand, so the same type also expresses the SV junctions of section 5.

### Coordinate lifting

`fragments.jl` owns `lift(frag, span)`, the only offset-to-reference conversion.
- `start` follows the prefixHashDB convention (`normalize_to_PAMseqEXT`, `db_helpers.jl:166`): the PAM-side end, i.e. `stop` for a forward `extends5` hit and `start` for reverse, flipped for `extends5 = false`. `end` is the other end.
- A coordinate that falls in a `RefSeg` maps through that segment, reversed when `isplus = false`.
- A coordinate that falls in a `LitSeg` maps to the `POS` of its insertion, and the allele is marked `ins` in the `variants` column. A site entirely inside a long insertion then has `start = end = POS`.
- The reported strand is the scan strand, flipped when the segment holding the PAM-side end has `isplus = false`.

### prefixHashScan integration

- **Variant work items.** The descriptors are cut into about 2 MiB items, appended as extra `PrefixHashScanChunkWork` entries after the genome items (`streaming.jl:356`). A worker builds the bytes of its item into its own buffer from FAI/2bit range reads plus the literals, then runs the existing kernel. Nothing loads a whole chromosome. The scan kernel needs no change.
- **Separator.** Fragments are separated by `distance + 1` copies of the byte `-`. It has IUPAC mask 0 (`verification.jl:3`) and raw code `0xff` (`kernel_common.jl:9`), so a window containing it is never a candidate, even with `ambig_max > 0`, and Myers scores it as a mismatch. An `N` run would match every guide base in Myers and send separator windows to the ambiguous slow path. With `-`, the records are about 30-50% smaller.
- **The anchor test runs before Myers**, and therefore before the early-stopping counter and before any row is written. Flank windows are plain reference sites and must not be reported twice or consume the variant cap.
- Fragment alleles are concrete bases. IUPAC appears in a fragment only from a real reference ambiguity code in a flank, or when the safety valve collapses SNPs deliberately. That removes vcfDB bug 3.
- A sorted vector of fragment start offsets per item maps each candidate to its fragment with one `searchsortedlast`.
- **Early stopping uses separate caps** for reference and variant hits, both defaulting to the same `early_stopping` values. `PrefixHashScanEarlyStopState` (`db_prefix_hash_scan.jl:274`) gets one `active_mask` and one counter set per record class, and a worker reads the mask of its item's class. "All guides retired" ends claiming only for that class, so variant items stay claimable after every guide fills its reference cap. The variant cap counts rows, i.e. walks, not loci: one locus with eight walks uses eight.
- **Engines.** `--vcf` requires a streaming backend. On the `:legacy` fallback (`hash_len` other than 16, motif span above 65) it is a startup error.

### Genotype modes

Section 7 states what each mode means for the output. The mechanics:

- `:all` (default): ignore GT and accept every walk. Correct for population VCFs, somatic and mosaic calls, and any VCF without genotypes.
- `:genotype`: accept a hit only if one sample could carry its sequence. The test runs in two places:
  - **at enumeration**, on the taken variants only: prune a walk when no sample carries all of its ALTs. Skipped variants are not tested here. The rule is monotone (a sample that carries a set carries every subset), so pruning never removes a walk that a smaller consistent walk needs.
  - **per hit**, after lifting, on every variant whose footprint lies inside the aligned span: a taken variant needs its ALT, a skipped one needs REF. This result fills the `samples` column, and a hit with no carrier is dropped.
  Testing skipped variants across the whole fragment instead would be a false negative: a sample that is `1/1` at a left-flank variant would reject every fragment anchored to its right, including sites that never touch that variant.
- **REF counts as an allele**, index `0`. A sample called `1/1` does not carry REF there, so a hit whose span covers that variant without taking it is rejected for that sample. Treating REF as always available would widen `:genotype` back toward `:all`.
- **Missing calls match anything.** `.` in a GT slot (`./.`, `./1`) matches every allele. Otherwise a single no-call would reject every walk across that variant for that sample, which is a false negative.
- **Phase is per sample and per record.** Two calls count as phased together only when both use `|` and, if `PS` is present, their `PS` values match. Otherwise the unphased rule applies to that pair.
- Phased: some haplotype slot of one sample carries every required allele. Unphased: every required allele appears somewhere in one sample's GT (cis assumed, the safe direction).
- **No ploidy option.** Both rules read the arity of each record as written, so haploid, diploid, mitochondrial and polyploid calls go through the same code.
- **Implementation.** Each allele gets a `BitVector` over haplotype slots (samples for the unphased rule). Consistency is a non-empty AND across the required alleles, and the same AND gives the carriers for the `samples` column. This keeps 2,504-sample cohorts cheap.

### Combinatorial safety valve

The cap is a cost trigger, not a limit on what gets searched. When an anchor's leaf count exceeds it:

1. Collapse the right-hand SNPs at one position into one IUPAC base covering REF and every ALT there, most-allele-first, until the count falls under the cap or `ambig_max` positions are spent. Each collapse removes those variants from the branching, and the IUPAC base is a superset of their alleles, so nothing is lost.
2. If the count is still over the cap, **emit every remaining walk anyway** and log the anchor with its count.

So the cap changes how a dense window is represented, never whether it is searched. Indels are never collapsed, because no IUPAC code expresses a length change.

The valve spends the motif's `ambig_max` budget, which defaults to 0 (`motif.jl:277`), so by default it never fires. Raising `ambig_max` also admits ambiguous reference windows in the genome scan. The valve is disabled in `:genotype` mode: a collapsed position has no single allele to test against a GT.

### Expected cost

| Input | Variants | Spacing | Variants per W window | Extra scan |
|---|---|---|---|---|
| Patient WGS | 4-5M | ~1 / 600 bp | 0-1 | ~5-10M fragments, 0.3-0.6 Gbp, roughly 10-20% of a genome scan |
| 1000G, common only | ~10M | ~1 / 300 bp | 0-2 | a few times the above |
| gnomAD v4 genomes | 700M+ | 1 / 4-8 bp in places | 4-8 | 16-256+ walks per anchor; many times the genome |

Numbers to confirm by measurement. A whole-genome gnomAD run stays the long-term target, but the supported workflow for population-scale input is a region-restricted run, with an upstream `bcftools view -i` frequency filter and the IUPAC valve as further levers. Section 8 defines the chr22 measurement.

## 5. Structural variants

In scope: `DEL`, `DUP`, `INV`, `BND`, and `INS` whose inserted sequence is in ALT or in `INFO/SEQ`.
Deferred with a warning and a count: symbolic `<INS>` and `<INS:ME:*>` without a sequence, `<CNV>`, and any record with a non-zero `CIPOS` or `CIEND`.

Key observation: an SV creates **new** sequence only at its breakpoint junctions and inside inserted sequence. The interior of a deletion or duplication is reference sequence the normal scan already covers on both strands, so it contributes no new sites. A copy number change changes a count, not a sequence.

Every supported SV therefore reduces to junction fragments in the same `Fragment` form as section 4. `W` is the window from section 4, and each side contributes `J = W - 1` bases, so any site spanning the junction fits. For symbolic records, `POS` is the padding base before the event (VCF 4.2), so the affected bases are `POS+1..END`. In the table, `Ref(c, a, b, ±)` is the inclusive range `a..b`, read reverse-complemented when `-`.

| SV | Junctions | Fragment segments | Anchor |
|---|---|---|---|
| `DEL`, POS=l, END=r | 1 | `Ref(c, l-J+1, l, +) + Ref(c, r+1, r+J, +)` | the 2 bases either side of the join |
| `DUP`, POS=s, END=e | 1 | `Ref(c, e-J+1, e, +) + Ref(c, s+1, s+J, +)` | the 2 bases either side of the join |
| `INV`, POS=s, END=e | 2 | left: `Ref(c, s-J+1, s, +) + Ref(c, e-J+1, e, -)`; right: `Ref(c, s+1, s+J, -) + Ref(c, e+1, e+J, +)` | the 2 bases either side of each join |
| `BND` | 1 per junction | `Ref(cA, ..., strandA) + [Lit(ins)] + Ref(cB, ..., strandB)`, from the ALT bracket notation | the join, plus any inserted bases |
| `INS` with sequence | 1 (counts as 2 joins) | `Ref(c, p-J+1, p, +) + Lit(seq) + Ref(c, p+1, p+J, +)` | the inserted span |

A deletion or insertion short enough to be written as a plain REF/ALT pair is already handled as an indel in section 4. The SV path takes over only for symbolic records, so there is one code path per representation, not per size.

### BND orientation

The four ALT forms in VCF 4.2 §5.4 map directly to segment strands. With the breakend at `cA:posA` and the mate at `cB:posB`:

| ALT | Meaning | Segments |
|---|---|---|
| `t[cB:posB[` | piece after `posB`, forward | `Ref(cA, ..posA, +) + Ref(cB, posB.., +)` |
| `t]cB:posB]` | piece before `posB`, reverse-complemented | `Ref(cA, ..posA, +) + Ref(cB, ..posB, -)` |
| `[cB:posB[t` | piece after `posB`, reverse-complemented, joined before `posA` | `Ref(cB, posB.., -) + Ref(cA, posA.., +)` |
| `]cB:posB]t` | piece before `posB`, forward, joined before `posA` | `Ref(cB, ..posB, +) + Ref(cA, posA.., +)` |

Bases between `t` and the bracket (`tAGT[cB:posB[`, `]cB:posB]AGTt`) are inserted at the junction and become a `LitSeg` between the two `RefSeg`s. Single breakends (`t.`, `.t`) have no mate sequence: they are skipped and counted.

Mate records (`MATEID`) describe the same junction from the other side, so emitting a fragment for both would double-report the junction. Deduplicate by the canonical junction key `(min, max)` of the two `(chrom, pos, side)` triples, and keep the record whose `ID` sorts first for the annotation.

### Small variants near a junction

The anchor rule from section 4 extends across a junction: a small variant within `W` of the breakpoint on either side joins the depth-first walk for that junction fragment. The only new part is that the walk-distance measurement crosses from one `RefSeg` to the next. Implement this after SNV, indel and plain junctions work, because it multiplies the test matrix. Until then, junction flanks are reference only, so a site that needs both a junction and a nearby small variant is missed. Section 10 lists this as a known gap.

### Reporting

A junction hit has two reference sides and no single `start`. `chromosome`, `start` and strand come from the segment holding the PAM-side end, by the lift rule in section 4; `end` comes from the segment holding the other end and may lie on another chromosome. Name both sides and name both sides in the `variants` column (`id|BND chrA:posA~chrB:posB`). Downstream `filter_overlapping` must not merge a junction hit with a reference hit at the same coordinate, because they are different sequences.
## 6. Output format

vcfDB wrote `guide,alignment_guide,alignment_reference,distance,chromosome,start,strand,variants`.

One merged table where reference hits leave `variants` empty. The new columns appear only when `--vcf` is given, so plain prefixHashScan output keeps its current schema. `filter` and `summarize` accept both schemas.

| Column | Notes |
|---|---|
| `guide` | unchanged |
| `alignment_guide`, `alignment_reference` | unchanged; the reference side shows the variant sequence |
| `distance` | unchanged |
| `chromosome`, `start` | unchanged, same position convention as prefixHashDB |
| `end` | new. The other end of the site by the lift rule (section 4). Indels change the site span, so `start` alone is not enough |
| `strand` | unchanged |
| `variants` | new format: only the alleles inside the aligned span, `;`-separated, each `id\|chrom:pos:REF>ALT` followed by any `--vcf_annotate` fields in the order given. Junctions use `id\|BND chrA:posA~chrB:posB` |
| `overlapping_variants` | new, reference rows only: IDs of VCF variants whose footprint lies inside the site, PAM included, `;`-separated. Empty on variant rows. Found by an interval lookup on the sorted variants, with no extra alignment |
| `samples` | `:genotype` mode only. Default: carrier count over samples with a call, e.g. `3/2504`. `--vcf_samples list` writes the comma-separated carrier names instead |

Rejected: a `ref_distance` column with the best reference-only distance at the same locus. It needs an extra alignment per variant hit and does not pay for itself. A variant that creates a PAM already shows up as a variant row with no reference row nearby, and `overlapping_variants` covers the opposite case: a reference row that lists variants but has no variant row next to it is a site the variants destroy or push past the search distance.

In `:genotype` mode reference rows are not filtered by genotype. A site where every sample is `1/1` for a variant inside it is still reported as a reference row, with that variant in `overlapping_variants`.

### Variant-aware summarize

`summarize_offtargets` (`src/utils.jl:481`) groups by guide and counts per distance. It currently counts every variant combination of one locus as a separate off-target, which inflates the counts in exactly the dense regions where the numbers matter.

Needed changes:
- use the locus definition of `filter_overlapping`: rows of one `(guide, chromosome, strand)` whose `start` values lie within `distance` of a kept row form one locus, which takes the minimum distance. `summarize --variant_aware` runs this grouping itself, so it gives the same answer on filtered and unfiltered input. **`end` stays out of the key.** A deletion inside a site shortens its reference span, so a variant row and the reference row at the same locus carry different `end` values. `start` is the PAM-side coordinate under the existing convention
- add `D0_var .. Dk_var`, counting loci whose lowest-distance row is a variant row, i.e. loci that the variants create or make closer
- add `D0_ref .. Dk_ref` for loci whose lowest-distance row is a reference row, and keep the existing `D0..Dk` as the union, so old output stays readable
- add `lost`, counting reference loci that list `overlapping_variants` and have no variant row: sites the variants may destroy
- in `:genotype` mode, add an optional per-sample or per-population breakdown behind a flag. Without it the table stays one row per guide
- CLI: `CHOPOFF summarize --variant_aware`

`filter_overlapping` already groups by locus in the way this needs; see below. The one change it does need is to not merge a junction hit with a reference hit at the same coordinate, since those are different sequences.

### Counts mode is reference-only

`prefixHashScan` has a second output mode that writes `guide,D0..Dk,complete` and skips traceback. **`--vcf` requires detail mode and errors in counts mode.**

The reason is not the anchor test, which works fine without traceback: the anchor and the candidate window are both ranges inside the fragment, so overlap is a range intersection. The reason is that counts mode has no rows to deduplicate. It increments per-distance counters during the scan and never materializes a locus. One locus with a reference site and eight variant walks would add nine to the count, so the number would measure walks rather than sites, and no later step could repair it.

The detail path already solves this with tools that exist:

```
CHOPOFF search ... prefixHashScan --vcf variants.vcf.gz   # detail rows
CHOPOFF filter --detail_file ... --distance 3             # collapse loci
CHOPOFF summarize --variant_aware                         # counts per guide
```

`filter_overlapping` (`src/utils.jl:425`) already does most of the collapsing this needs. It sorts by `(guide, chromosome, strand, distance, start)` and drops any row whose `start` falls within `distance` of a row already kept, so a reference row and the variant rows at the same locus collapse to the single lowest-distance row. One change and one consequence:

- **A variant row is dropped when the reference row is as good or better.** That is correct for safety: if the site is already at distance 2 without any variant, the variant adds no risk there. A variant that *improves* the match, including one that creates a PAM, has the lower distance and is the row that survives.
- **On a tie, keep the reference row.** The site exists without any variant, so that is the honest representative. This needs a code change: the current sort key lets a variant row with a smaller `start` win a tie. Add `!isempty(variants)` to the sort key after `distance` and before `start`. The kept reference row still carries `overlapping_variants`, so the variant at that locus stays visible.

One scaling note, not a blocker: `filter_overlapping` builds a `Set` of `2 * distance + 1` integers per kept row and a joined string key per row. A variant run makes the detail table larger, so this becomes the slowest step in the pipeline well before the search does. Worth measuring during the chr22 run, and worth replacing with a sorted sweep if it bites.

## 7. Frequencies, samples, and what a row means

CHOPOFF does no frequency filtering and no sample selection. It consumes the VCF as given, searches every allele in it, and copies the fields the user names into the output. Selecting which variants are in play is the user's step, done before the run with `bcftools`.

### Why filtering belongs upstream

- There is no single AF field. gnomAD ships `AF`, `AF_afr`, `AF_nfe`, `AF_grpmax` and more; somatic callers put VAF in `FORMAT/AF` or compute it from `FORMAT/AD`; many VCFs have no frequency field at all. A `--min_af` flag would have to guess which field the user means, and guessing wrong silently drops alleles. That is the one failure class this feature exists to prevent.
- `bcftools view -i` already expresses every filter anyone would want, including ones CHOPOFF could not: per-population thresholds, AC and AN combinations, `FORMAT/DP` floors, sample subsets, region masks.
- A frequency says nothing about whether the genome in front of you carries the allele. Filtering on AF discards exactly the private and low-frequency alleles that a personalized search exists to find.

So `min_af` and `pass_only` are dropped from the design. They were two flags, two code paths and a class of silent false negatives, in exchange for one `bcftools` call the user can write better.

### What the user decides instead

The real knob is not frequency, it is which VCF and which mode. Those two choices fix the meaning of every output row.

| Input | Mode | A variant row means |
|---|---|---|
| one patient, phased | `:genotype` | this patient has this site on this haplotype |
| one patient, unphased | `:genotype` | this patient has this site, if the alleles are in cis |
| one patient, somatic/mosaic | `:all` | some cell population may have this site |
| multi-sample cohort | `:genotype` | at least one listed sample has this site; the `samples` column names them |
| multi-sample cohort | `:all` | some combination of cohort alleles gives this site, possibly a chimera no sample carries |
| gnomAD or other AF-only | `:all` | some human genome may have this site. Not a claim about any individual |

The chimeric case is worth stating plainly in the user docs: `:all` on a merged multi-sample VCF is a superset, so it is safe, but it reports haplotypes that may exist in nobody. A user who wants per-individual truth runs `:genotype`, or splits the cohort with `bcftools view -s SAMPLE` and runs each sample. CHOPOFF does not merge or normalize cohorts, which is a well-served job for other tools.

### What each mode does, mechanically

Both modes start from the same walks (section 4). `:genotype` adds two filters: it prunes walks whose taken ALTs no sample carries, and it checks each hit against the alleles inside its span.

`:all` accepts every combination.
- It never looks at `FORMAT`, so it works on any VCF, including one with no sample columns.
- The result is the union over every genome the VCF could describe. It is a superset of what any single genome contains.
- Nothing is lost. Combinations that no real haplotype carries are reported as extra rows, which is the safe direction.

`:genotype` accepts a combination only if one sample could carry all of its alleles together.
- Phased, same `PS` block: the alleles must sit on one haplotype of one sample.
- Phased, different `PS` blocks, or unphased: the alleles must all appear in that sample's genotype, in any arrangement. Phase is unknown, so cis is assumed, which stays on the safe side.
- Arity is read per record, so haploid, diploid and polyploid calls need no configuration.
- Every `:genotype` site is also an `:all` site. Rows can differ only where `:all` used the IUPAC valve, which `:genotype` never uses. The result is smaller and more specific, and it is only correct when a sample column is a real individual.

The rule to state in the user docs: `:all` answers "can any genome have this site?", `:genotype` answers "can one of these named individuals have this site?". Choose `:genotype` only when the sample columns are individuals you care about.

### gnomAD with `:genotype`

Two outcomes, depending on which file the user actually has.

**A real gnomAD sites-only VCF: startup error.** gnomAD releases carry the eight mandatory columns and no `FORMAT` or sample columns. There is no genotype to test, so `--vcf_mode genotype` fails immediately:

```
ERROR: --vcf_mode genotype requires per-sample genotypes, but
       gnomad.chr22.vcf.gz has 0 sample columns and no FORMAT field.
       gnomAD is a frequency resource, not a set of individuals. Use
       --vcf_mode all, and restrict to a population upstream if needed:
         bcftools view -i 'INFO/AF_afr>0' gnomad.chr22.vcf.gz
```

This is the behaviour we want: loud, immediate, and it names the fix.

**A gnomAD file converted to pseudo-samples: it runs, and the output is wrong.** This is the real hazard. CRISPRme's gnomAD converter turns each population into a fake sample, so the file gains GT columns named `afr`, `nfe`, `eas` and so on. CHOPOFF cannot tell those from individuals, so `:genotype` succeeds and produces confident nonsense:

- A combination is kept when all of its alleles appear in the same population. Co-occurrence in a population of 76k people says nothing about co-occurrence in one person, so kept rows overstate what they show.
- Worse, combinations whose alleles sit in different populations are **dropped**. An admixed individual can carry exactly those, so this loses real haplotypes. The mode becomes restrictive without becoming correct, which is the one direction this design refuses to take.
- The `samples` column would then name populations, reading as if `afr` were a person.

So: on pseudo-sample input, `:genotype` is a false-negative generator. `:all` on the same file is correct and is what the data supports.

Detecting pseudo-samples reliably is not possible from the VCF alone. The cheap partial guard is to warn in `:genotype` mode when no sample carries any per-sample depth-like field (`DP`, `AD`, `GQ`), since real individual calls almost always carry one and converter output carries none. That is a heuristic warning, never an error. Deferred unless it proves needed, with the documentation carrying the weight in the meantime.

### The on-target site needs no special handling

CHOPOFF searches for every sequence within the given distance of the guide. The intended target is one of those sequences, at distance 0, so it is already in the output. Variants inside it are covered by the same walks as any other locus, which means the variant-aware search is a superset of a variant-aware on-target analysis. No on-target mode, no second code path, no separate report.

What the user reads off the **unfiltered** table, by grouping rows at the on-target locus (`filter` keeps one row per locus, so it removes the distinction between the first three rows):

| Rows at the on-target locus | Meaning |
|---|---|
| reference row at distance 0, empty `overlapping_variants` | the reference target; no variant lies in it |
| reference row at distance 0 plus a variant row at distance 0 | a variant sits in the site but leaves the match intact |
| reference row at distance 0 plus a variant row at distance > 0 | a variant changes the match, and the distance says by how many edits |
| reference row at distance 0 that lists `overlapping_variants`, and **no** variant row | a variant destroys the PAM or pushes the site past the search distance, so a carrier may not be cut at all |

CHOPOFF reports edit distance, not predicted activity. It does not score how much a variant reduces cutting, which is what CRISPR-HAWK adds on top. The claim here is narrower and checkable: this sequence is or is not within the searched distance of the guide.

The last line still depends on a missing variant row, but `overlapping_variants` puts the evidence on a row the user can see, and `summarize --variant_aware` counts these loci in `lost`. In `:genotype` mode, also check `samples` on the variant rows: a heterozygous carrier keeps the reference site on the other haplotype.

### Annotation passthrough

Instead of filtering, CHOPOFF copies fields verbatim:

```
--vcf_annotate INFO/AF,INFO/AF_nfe,FORMAT/AF
```

- Default is the record `ID` alone, or `CHROM:POS:REF>ALT` when `ID` is `.`.
- Named `INFO` keys are appended to each allele in the `variants` column, in the order given.
- Named `FORMAT` keys resolve per sample, so they only make sense together with `:genotype` mode or a single-sample VCF. With many samples the field grows, which is the user's choice to make, not a default.
- A named key that is absent from the header is an error at startup, not a silent empty column.

This keeps one rule: whatever the user wants to rank or threshold on afterwards, they name it and it appears in the table. Ranking and thresholding then happen in `summarize`, `filter`, or any downstream script.

## 8. Test plan

Tests run on the existing example genome (`test/sample_data/genome/semirandom.fa`, 8 chromosomes of 108326 bp) with an upgraded VCF. Everything stays inside the existing `test/src/` auto-include mechanism.

### Test assets

1. **Replace `test/sample_data/artificial.vcf`.** It is unusable as a correctness reference because it lists REF inside ALT (`A -> A,G,T`), which hides vcfDB bug 1, and it is not normalized. Keep it under a new name for the vcfDB regression comparison and write a new file set.
2. **New files** under `test/sample_data/vcf/`:

| File | Contents | Purpose |
|---|---|---|
| `small_norm.vcf` | normalized SNVs and indels, one ALT per line, on `semirandom1` | the main correctness file |
| `small_phased.vcf` | the same sites with 4 phased samples and 2 `PS` blocks | `:genotype` phased, including the block boundary |
| `small_unphased.vcf` | the same sites with unphased GT | `:genotype` unphased |
| `mixed_arity.vcf` | haploid, diploid and hexaploid GT in one file | arity read per record, no option |
| `population.vcf` | `INFO/AF` only, no samples or FORMAT column | `:all` plus AF annotation; `:genotype` must error |
| `sv.vcf` | one each of DEL, DUP, INV, BND pair, INS with sequence | section 5 |
| `malformed.vcf` | multi-allelic lines, unsorted POS, REF mismatch, `ALT=.`, `ALT=*`, `<CNV>`, `<INS>`, non-zero CIPOS, unknown chromosome, non-left-aligned indel | the section 3 contract, one case per expected warning or error |
| `wrong_build.vcf` | every REF mismatching, and `chr`-prefixed names against a plain-numbered reference | the two hard-error paths |

Each site in `small_norm.vcf` is placed deliberately against a known PAM in `semirandom1`, not at random, so the expected hit set is small and reviewable. Generate them with a committed script (`test/sample_data/vcf/make_test_vcfs.jl`) so the placement logic is inspectable and regenerable.

### Correctness cases

Every case below needs an expected hit set from the oracle, not a golden file produced by the implementation.

| Case | Why |
|---|---|
| two SNPs at every spacing 1..`lp` | vcfDB bug 2 (window too short) |
| a SNP and an indel in one window | the combination CRISPRme cannot do |
| REF+ALT mixed combinations of two sites | vcfDB bug 1 |
| a SNP inside a deletion | vcfDB bug 5 (mutually exclusive alleles) |
| two overlapping indels | conflicting variants never share a walk |
| deletion 10-20, SNP at 12, deletion 18-25 | SNP plus second deletion must combine (conflict is not transitive) |
| two insertions at the same gap | exclusive; each appears alone |
| a deletion inside a site | the anchor is the join; the site must be found |
| a site ending inside an insertion | `start` lifts to the insertion `POS`, allele marked `ins` |
| an insertion longer than the site | vcfDB bug 4 |
| a multi-ALT site, pre-split and not | the contract path |
| a SNP at `p` and an insertion after `p` | both must combine; the gap-edge rule in section 4 |
| a window overlapping a `-` separator, and a flank window overlapping no anchor | must yield no row and must not consume the variant cap |
| a site ending exactly `W - 1` past its anchor | the completeness invariant; nothing cut off at a fragment edge |
| a sample called `1/1` at a variant | `:genotype` must reject a hit whose span covers it as REF |
| a sample `1/1` at a left-flank variant, site not touching it | `:genotype` must still report the site (per-hit test, not per fragment) |
| a no-call (`./.`, `./1`) at a variant inside a site | matches any allele; the site must be reported |
| a phased file with one `/` call | that pair falls back to the unphased rule |
| a reference row with a PAM-destroying variant | `overlapping_variants` lists it; `summarize` counts the locus in `lost` |
| a guide that retires on the reference | variant work items are still claimed (separate active masks) |
| `--vcf` with a `:legacy`-only configuration | startup error |
| a deletion inside a site, through `summarize` | the variant and reference rows must group to one locus despite different `end` |
| a VCF with `chr1` naming against a `1` reference | must error, not report an empty variant set |
| a VCF called against the wrong build | above the 1% REF-mismatch rate, must error |
| `--vcf` with counts mode | must error, not silently count walks |
| a locus with a reference row and several variant rows | `filter` must collapse them to the lowest-distance row |
| a variant that improves the match | the variant row must survive the filter, not the reference row |
| a haploid GT (`1`) on a sex chromosome | per-record GT arity |
| a variant within `lp` of a chromosome edge | vcfDB bug 10 |
| a variant that creates a PAM where the reference has none | the main clinical motivation |
| a variant that destroys a PAM | the hit must disappear from the variant walk but stay as a reference hit |
| reverse-strand sites | vcfDB bug 7 |
| `extends5 = false` (Cas12a) | vcfDB bug 8; mirrors `test/src/db_extends5_false.jl` |
| distance 0 through 4 | the scan envelope |
| bulges at the variant boundary | insertion and deletion interaction with the alignment |
| `:genotype` phased, inside and across a `PS` block | the section 3 trap |
| `:genotype` unphased, haploid, diploid and hexaploid GT | arity-agnostic consistency rule |
| a guide that fills its reference cap | its variant sites must still be reported (separate caps) |
| the on-target site with a variant in it | must appear as a variant row at low distance |
| the IUPAC safety valve, forced by a low cap | the fallback must be a superset, never a subset |
| `:genotype` on a sites-only VCF | must be a startup error, not an empty result |
| `:all` and `:genotype` on the same cohort file | every `:genotype` site must be covered by an `:all` row (valve off: a subset) |
| each SV type from `sv.vcf` | section 5 |
| `--region` with `--vcf`, a variant just outside the region inside a site that overlaps it | the variant is used and the row reported |
| `--region` cutting a BND pair | the junction is kept when either breakend is near the region |
| a BND mate pair | junction deduplication |

### Oracle

The reference implementation for the tests, kept in the test suite and never in `src/`:

```
for every combination of alleles in the region:
    build the full alternative sequence of that haplotype
    write it as a one-chromosome FASTA
    run search_prefixHashScan with query_variant=:bruteforce on it
    lift each hit back to reference coordinates with its own offset map
union the hits
```

The bruteforce variant runs on the legacy engine and skips the prefix filter, and the oracle lifts coordinates with its own offset map, so it shares no enumeration, lifting or scan code with the engine under test.

In `:all` mode the new engine must return this set. A missing hit fails the test. An extra hit is reported, and it is allowed only when the IUPAC valve was used or when the candidate-window anchor test (section 4) let a non-anchor alignment through; the test classifies which.

The oracle is exponential, so it runs only on small windows of `semirandom1`. That is why the test VCF is hand-placed rather than dense.

### Scale tests

Not part of `runtests.jl`. These live under the ignored `test/local_human/` tree, following the existing benchmark layout.

- a region run: `--region chr22` on GRCh38 with the gnomAD chr22 VCF, measuring fragment count, extra bases scanned, and wall time against the same region without `--vcf`
- a patient-scale run: one 1000G sample, whole genome, `:genotype` phased
- record both in `prefixHashScan.md` style, with medians over alternating runs, because the host is shared

### Regression against vcfDB

Run the old `artificial.vcf` with REF removed from the ALT lists through both engines and record which sites vcfDB misses. This is a one-time document in this file, not a permanent test, since vcfDB is being retired.

### Separate test for bug 12

A site whose first base is ambiguous, read back through `search_prefixHashDB`, to confirm or rule out the `sum(is_ambig[1:bv_start])` off-by-one in `src/db_prefix_hash.jl:445`. This is independent of the VCF work.
## 9. Code structure

The feature is large, so it is split by responsibility and each part is testable alone. Nothing here changes the scan kernel.

```
src/variants/
  vcf_read.jl       parse records, validate the section 3 contract, emit Variant
  variants_index.jl footprints, anchor spans, walk order, conflict test, interval lookup
  genotypes.jl      GT, phasing, PS blocks -> per-allele carrier bitsets
  fragments.jl      Fragment, RefSeg, LitSeg, materialize, lift coordinates
  enumerate.jl      leftmost-ALT anchor, depth-first walk, combination cap, IUPAC valve
  sv.jl             SVTYPE/END/SVLEN/BND -> junction fragments
  variant_record.jl descriptor work items, record bytes, `-` separators, offset table, anchor test
src/variants.jl     the module entry, re-exported from src/CHOPOFF.jl
```

Boundaries that keep this clean:
- `vcf_read.jl` is the only file that knows the VCF format. Everything downstream works on `Variant`. It reads plain and bgzipped VCF with the existing `VariantCallFormat.VCF.Reader` over a `BGZFStreams` stream; both packages are already in the Manifest.
- `fragments.jl` owns coordinate lifting. No other file converts fragment offsets to reference positions, which is where vcfDB accumulated four inconsistent formulas (bug 8).
- `genotypes.jl` exposes `carriers(alleles, mode)`, the AND of per-allele carrier bitsets. Enumeration prunes on a non-empty result for the taken alleles; the per-hit test uses the span's alleles. `:all` mode never calls it, so the hot path has no genotype branch.
- `sv.jl` only produces `Fragment`s. The enumerator and the search path cannot tell an SV junction from an indel.
- `variant_record.jl` is the only point of contact with `search_prefixHashScan`. It owns the variant work items, record building, the offset table, and the anchor test that gates every hit. The scheduler change it needs is the per-class early-stop state (section 4).

Suggested order of work, each step independently testable:

1. `vcf_read.jl` plus the contract tests on `malformed.vcf`
2. `variants_index.jl` plus `fragments.jl`, tested by materializing and lifting without any search
3. `enumerate.jl` in `:all` mode, SNV and indel only, validated against the oracle
4. `variant_record.jl` and the prefixHashScan hookup, end to end on `small_norm.vcf`
5. output columns and `summarize --variant_aware`
6. `genotypes.jl` with the phased, unphased, `PS` and mixed-arity cases
7. `sv.jl` for DEL, DUP, INV, INS, then BND with mate deduplication
8. small variants near a junction
9. the scale runs on chr22

Steps 1 through 5 are a complete, useful feature on their own: a correct `:all` mode search over SNVs and indels. That is already more than vcfDB does.

### CLI

```
CHOPOFF search --database ... --guides ... --output ... prefixHashScan \
    --vcf variants.vcf.gz \
    --vcf_mode all|genotype \
    --region chr22 \
    --vcf_annotate INFO/AF,FORMAT/AF \
    --vcf_samples count|list
```

- `--vcf_mode genotype` without GT in the file is an error, not a silent downgrade.
- `--region` is a general prefixHashScan option (`prefixHashScan.md`, Planned: region restriction). It limits the reference scan and the variants together, so reference and variant rows cover the same loci. With `--vcf`:
  - Variants are kept when their footprint lies within `W` of a region, because a site overlapping the region can use a variant up to `W` bases outside it. No tabix reader exists in the dependencies, so this is a streaming filter over the whole file.
  - Variant rows use the same rule as reference rows: the candidate window, lifted to the reference, overlaps a region. Inside a `LitSeg` the lift gives the insertion `POS` (section 4).
  - A junction is kept when either breakend lies within `W` of a region, and its rows when any `RefSeg` part of the candidate window overlaps a region. One BND record describes its junction fully, so a mate record outside the region is not needed.
  - The section 3 contract checks apply to the kept records. The "no chromosome matches" error applies to the region's chromosomes.
- `--vcf_annotate` names fields to copy into the output (section 7). An unknown key is an error at startup.
- `--vcf_samples` (default `count`) sets the `samples` column format in `:genotype` mode.
- There is no frequency or FILTER flag. That selection happens upstream in `bcftools`.

### Julia API

```julia
search_prefixHashScan(dbi, guides, output; ...,      # region: see prefixHashScan.md
    vcf::Union{Nothing, String} = nothing,
    vcf_mode::Symbol = :all,               # :all or :genotype
    vcf_annotate::Vector{String} = String[],
    vcf_samples::Symbol = :count)          # :count or :list
```

`PrefixHashScanStats` gains `variants_read`, `variants_skipped` (per reason), `fragments`, `fragment_bases`, `valve_anchors`, `variant_candidates`, and `anchor_rejected`. With `verbose=true`, a `prefixHashScan variants` `@info` record after enumeration reports the same counts, plus the descriptor memory and the enumeration time. This record is the "run summary" that section 3 refers to.
## 10. Decisions and deferred work

Decided:
- Only `search_prefixHashScan` gets VCF support. `vcfDB` is retired, not fixed.
- No frequency filtering, FILTER handling or sample merging at all. CHOPOFF searches the VCF as given and copies the fields the user names (section 7). `min_af` and `pass_only` are dropped from the design.
- Early stopping uses separate caps for reference and variant hits, both defaulting to the same `early_stopping` values. A guide that saturates on the genome still reports its variant sites. This must be stated in the user docs, since the totals no longer add up to one shared budget.
- No ploidy option. The consistency rule reads GT arity per record, so haploid, diploid and polyploid calls all work unchanged.
- No on-target mode. The intended target is a distance-0 row like any other, so the variant search already covers it (section 7).
- No row budget. A run produces as many rows as the caps allow and takes as long as it takes.
- `--vcf` works in detail mode only. Counts mode would count walks rather than loci and has no row to deduplicate, so it errors instead. The per-locus counts come from the existing `search` then `filter` then `summarize` pipeline (section 6).
- Walks that reproduce reference sequence are kept and annotated. They are rare, and supplying a clean normalized VCF is the user's job.
- Users normalize their own VCFs. CHOPOFF validates, repairs only multi-ALT splitting and IUPAC expansion, and reports (section 3).
- Variants are independent on/off choices with a pairwise conflict test, not merged groups (section 4).
- The anchor test uses the candidate window and runs before Myers. Duplicate rows are allowed; missed rows are not.
- Fragments are descriptors built once per run and scanned as extra streaming work items, separated by `distance + 1` `-` bytes. `--vcf` on a `:legacy`-only configuration is a startup error.
- Early stopping keeps one active mask per record class.
- `:genotype` prunes on the taken alleles during enumeration, then checks each hit on its span. No-calls match any allele. The IUPAC valve is off in `:genotype` mode.
- Reference rows carry `overlapping_variants`. `filter` keeps the reference row on a distance tie. `summarize --variant_aware` uses the `filter` locus definition and reports `lost`.
- `samples` defaults to a carrier count; `--vcf_samples list` writes names.
- `--region` is a general prefixHashScan option that limits the reference scan and the variants together. A site is in a region when its candidate window overlaps it.
- `summarize` gets a variant-aware mode.
- `:genotype` on a VCF without sample columns is a startup error that names the fix. On pseudo-sample input (gnomAD converted per population) it cannot be detected reliably, so the documentation carries that warning (section 7).
- Correctness tests run on the example genome with hand-placed variants and an oracle. Scale runs use one chromosome or region of GRCh38.
- SVs in scope: DEL, DUP, INV, BND, and INS with a sequence.

Deferred, each with a warning and a count in the run summary:
- symbolic `<INS>`, `<INS:ME:*>` without a sequence. A mobile-element consensus approximation is possible later but would report sites that may not exist, so it needs an explicit opt-in flag.
- imprecise breakpoints (`CIPOS`, `CIEND`). Enumerating one fragment per offset in the interval is correct but the cost scales with the interval width, so it needs a cap design first.
- `<CNV>`: no new sequence, so nothing to search. It may still matter for on-target dosage, which is out of scope here.
- small variants within `W` of an SV junction (section 5, step 8). Until then, a site that needs both is a known false negative.
- an indexed (tabix) reader for `--region`. Add one only if the streaming filter is too slow on real files.
- whole-genome gnomAD: the design supports it, but the first release documents a region-restricted workflow and the measured cost rather than promising a genome-wide run.

- Per-sample VAF versus population AF in the output is resolved by `--vcf_annotate`: the user names the fields, so neither is privileged and the default stays short.

Documentation owes the user these, because none is visible from the output alone:
- the two independent caps, and that a saturated reference cap does not suppress variant rows
- how to read the on-target locus from the unfiltered table, in particular that a reference row with `overlapping_variants` and no variant row means a destroyed PAM (section 7)
- that the variant cap counts walks, not loci
- that the IUPAC valve needs `ambig_max > 0` and that this also widens the genome scan
- that `:genotype` does not filter reference rows by genotype
- that `:all` on a merged cohort reports haplotypes no individual carries, and that `:genotype` on pseudo-sample gnomAD loses real ones (section 7)
- that every frequency and sample decision happens upstream in `bcftools` (section 7)
- that non-diploid input needs no configuration
- that reported numbers are edit distances, not activity predictions
- that a variant run needs `filter` before `summarize`, and that `filter` keeps the lowest-distance row per locus, so a variant row disappears when the reference row is as good; the reference row keeps the variant in `overlapping_variants`


Open:
- None blocking implementation.
