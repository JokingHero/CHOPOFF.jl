# Record helpers shared by test/local_human/vcf_realistic.jl, the real-data
# fixtures (test/sample_data/vcf/real/) and their test. Independent of the
# production VCF reader. Needs test/helpers/vcf.jl.
using Random

const IUPAC = Dict('A'=>"A",'C'=>"C",'G'=>"G",'T'=>"T",'R'=>"AG",'Y'=>"CT",'S'=>"CG",'W'=>"AT",'K'=>"GT",
    'M'=>"AC",'B'=>"CGT",'D'=>"AGT",'H'=>"ACT",'V'=>"ACG",'N'=>"ACGT")

# Symbolic ALT as searched: "<TYPE>/END" for precise DEL*, DUP, DUP:TANDEM and
# INV (END from INFO/END, else POS + |SVLEN|), "" for skipped kinds, `nothing`
# for sequence ALTs. Kept in its original case, like the production allele key.
function symbolic_alt(alt, info, pos, alt_idx)
    any(c -> c in "<[].", alt) || return nothing
    type = uppercase(alt)
    supported = startswith(type, "<DEL>") || startswith(type, "<DEL:") || type in ("<DUP>", "<DUP:TANDEM>", "<INV>")
    items = split(info, ';')
    (supported && !("IMPRECISE" in items)) || return ""
    value(key) = (i = findfirst(x -> startswith(x, key * "="), items); i === nothing ? nothing : items[i][(length(key) + 2):end])
    stop = value("END") !== nothing ? parse(Int, value("END")) :
        (v = split(value("SVLEN"), ','); pos + abs(parse(Int, v[length(v) == 1 ? 1 : alt_idx])))
    return "$alt/$stop"
end

# Reference footprint [first, stop) and sequence of a searched allele.
function footprint(ref, alt, pos, seq = nothing)
    startswith(alt, '<') || return trim_footprint(ref, alt, pos)
    stop = parse(Int, split(alt, '/')[end])
    type = uppercase(alt)
    startswith(type, "<DEL") && return pos + 1, stop + 1, ""
    copy = seq === nothing ? "" : seq[(pos + 1):stop]
    startswith(type, "<DUP") && return pos + 1, pos + 1, copy
    return pos + 1, stop + 1, seq === nothing ? "" : vtest_rc(copy)
end

# "" marks ALTs that are not searched ("*", skipped SV kinds), so ALT numbers stay aligned.
searched_alt(a, info, pos, i) = a in (".", "*", "<NON_REF>", "<*>") ? "" :
    (sv = symbolic_alt(a, info, pos, i); sv === nothing ? uppercase(a) : sv)

function trim_footprint(ref, alt, pos)
    p = 0
    while p < min(length(ref), length(alt)) && ref[p + 1] == alt[p + 1] && ref[p + 1] in "ACGT"
        p += 1
    end
    s = 0
    while s < min(length(ref), length(alt)) - p && ref[end - s] == alt[end - s] && ref[end - s] in "ACGT"
        s += 1
    end
    p == length(ref) == length(alt) && ((p, s) = (0, 0))
    return pos + p, pos + length(ref) - s, alt[(p + 1):(end - s)]
end

function concrete_alleles(alt)
    out = [""]
    for c in alt
        out = [o * b for o in out for b in IUPAC[c]]
    end
    return out
end

haplotype_count(groups) = prod(1 + length(g) for g in groups; init = 1)

# Oracle edits of literal records, in coordinates relative to window start `ws`.
function window_edits(seq, records, ws, we)
    edits = NamedTuple[]
    for (k, r) in enumerate(records)
        for a in r.alts, c in (isempty(a) ? String[] : concrete_alleles(a))
            first, stop, alt = trim_footprint(r.ref, c, r.pos)
            push!(edits, (; first = first - ws + 1, stop = stop - ws + 1, alt, id = r.id, group = "g$k",
                pos = r.pos - ws + 1, ref = r.ref, vcf_alt = a, concrete = c))
        end
    end
    return edits
end

# Up to `n` guides from sites touching an edit (or any reference site), with up to `d` substitutions.
function guides_for(rng, w, edits, d; n = 3, haps = vtest_haplotypes(w, edits; chrom = "w"))
    guides = String[]
    for (hap, _, events, alternative) in shuffle(rng, haps)
        length(guides) >= n && break
        sites = Tuple{Int, Bool}[]
        for p in 1:(length(hap) - 22)
            hap[p + 21] == 'G' && hap[p + 22] == 'G' && push!(sites, (p, false))
            hap[p] == 'C' && hap[p + 1] == 'C' && push!(sites, (p, true))
        end
        near = filter(sites) do (p, rev)
            any(ev -> max(ev[1], p) <= min(ev[2], p + 22) || (ev[1] > ev[2] && p < ev[1] <= p + 22), events) ||
                !alternative
        end
        isempty(near) && continue
        p, rev = rand(rng, near)
        g = collect(rev ? vtest_rc(hap[p:(p + 22)])[1:20] : hap[p:(p + 19)])
        for _ in 1:rand(rng, 0:d)
            i = rand(rng, 1:20); g[i] = rand(rng, collect("ACGT"))
        end
        any(c -> !(c in "ACGT"), g) && continue
        push!(guides, String(g))
    end
    return unique(guides)
end

# GT/PS of one sample column: (alleles, label) with ALT numbers (0 REF, -1 missing)
# and label :chrom (`|` or haploid, no PS), :unphased or a PS integer.
function vreal_call(format, value)
    keys, v = split(format, ':'), split(value, ':')
    gi, pi = findfirst(==("GT"), keys), findfirst(==("PS"), keys)
    gt = v[gi]
    alleles = [a == "." ? -1 : parse(Int, a) for a in split(gt, r"[|/]") if !isempty(a)]
    ps = pi === nothing || length(v) < pi ? "." : v[pi]
    label = occursin('/', gt) && length(alleles) > 1 ? :unphased : ps == "." ? :chrom : parse(Int, ps)
    return alleles, label
end

# Records of a plain-text VCF; `alts` keep their ALT numbers ("" when not searched).
function vreal_records(path; contig = nothing)
    names = String[]
    records = NamedTuple[]
    for line in eachline(path)
        startswith(line, "##") && continue
        f = split(line, '\t')
        if startswith(line, '#')
            names = String.(f[10:end])
            continue
        end
        contig === nothing || f[1] == contig || continue
        pos = parse(Int, f[2])
        calls = Dict(name => vreal_call(f[9], f[9 + k]) for (k, name) in enumerate(names))
        push!(records, (chrom = String(f[1]), pos, id = String(f[3]), ref = uppercase(String(f[4])),
            alts = [searched_alt(a, f[8], pos, i) for (i, a) in enumerate(split(f[5], ','))],
            calls, fields = String.(f)))
    end
    return names, records
end

# Literal REF and ALT of a searched symbolic DEL, DUP or INV at `pos` (1-based in `seq`).
function vreal_literal(seq, pos, alt)
    first, stop, inserted = footprint("", alt, pos, seq)
    return seq[pos:(stop - 1)], seq[pos:pos] * inserted
end
