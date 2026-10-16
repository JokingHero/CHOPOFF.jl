# Symbolic DEL, DUP and INV are searched as their literal equivalents. The literal
# path is checked against the oracle; the symbolic run must give the same rows
# once allele keys are renamed.

# Literal edit and symbolic record for an SV over reference bases pos+1:stop.
function vtest_sv(r, kind, pos, stop, id; alt = "<$kind>", info = "END=$stop")
    edit = kind == "DEL" ? vtest_edit(r, pos + 1, stop + 1, "", id) :
        kind == "DUP" ? vtest_edit(r, pos + 1, pos + 1, r[(pos + 1):stop], id) :
        vtest_edit(r, pos + 1, stop + 1, vtest_rc(r[(pos + 1):stop]), id)
    record = ("chr1", pos, id, r[pos:pos], alt, ".", "PASS", info)
    return edit, record, "chr1:$pos:$(r[pos])>$alt/$stop"
end

vtest_sv_rows(path) = Set(vtest_signature(row) for row in eachrow(DataFrame(CSV.File(path))))

# Runs the oracle-checked literal search and the symbolic search on one case.
function vtest_sv_compare(dir, r, svs, guides, motif, distance; extra = Tuple[], chunk = 64, header = String[])
    edits = [first(s) for s in svs]
    literal = vtest_compare(dir, r, edits, guides, motif, distance; chunk, case_id = "sv_literal")
    genome = joinpath(dir, "reference.fa")
    records = sort!(vcat([s[2] for s in svs], extra); by = x -> x[2])
    vcf = vtest_vcf(joinpath(dir, "symbolic.vcf"), records; header)
    out = joinpath(dir, "symbolic.csv")
    stats = PrefixHashScanStats()
    CHOPOFF.search_prefixHashScan(guides, genome, motif, out; distance, vcf, stream_chunk_bases = chunk,
        early_stopping = fill(typemax(Int), distance + 1), vcf_max_alleles = vtest_exact_cap(edits), stats)
    renamed = joinpath(dir, "literal_renamed.csv")
    text = read(joinpath(dir, "actual.csv"), String)
    for (edit, _, key) in svs
        text = replace(text, vtest_key(edit, "chr1") => key)
    end
    write(renamed, text)
    @test vtest_sv_rows(out) == vtest_sv_rows(renamed)
    return DataFrame(CSV.File(out)), stats
end

vtest_sv_random(rng, n) = String(rand(rng, collect("ACGT"), n))

# Random case where the SV joins two halves of a Cas9 site at a random split:
# DEL removes the bases between them, DUP repeats a region ending with the first
# half and starting with the second, INV turns the second half around. The
# guide gets up to `d` substitutions; up to three SNPs are added.
function vtest_sv_random_case(rng, i)
    g = vtest_sv_random(rng, 20)
    site = g * "AGG"
    k = rand(rng, 1:20)
    kind = ("DEL", "DUP", "INV")[mod1(i, 3)]
    flank() = vtest_sv_random(rng, rand(rng, 30:50))
    middle = vtest_sv_random(rng, rand(rng, 1:60))
    local r, pos, stop
    while true
        left = flank()
        if kind == "DEL"
            r = left * site[1:k] * middle * site[(k + 1):end] * flank()
            pos = length(left) + k
            stop = pos + length(middle)
        elseif kind == "DUP"
            region = site[(k + 1):end] * middle * site[1:k]
            r = left * region * flank()
            pos, stop = length(left), length(left) + length(region)
        else
            region = middle * vtest_rc(site[(k + 1):end])
            r = left * site[1:k] * region * flank()
            pos = length(left) + k
            stop = pos + length(region)
        end
        # A literal inversion with complementary ends would be trimmed to a smaller footprint.
        kind == "INV" && r[pos + 1] == only(vtest_rc(r[stop:stop])) && (middle = vtest_sv_random(rng, rand(rng, 1:60)); continue)
        break
    end
    d = rand(rng, 0:3)
    guide = collect(g)
    for _ in 1:rand(rng, 0:d)
        guide[rand(rng, 1:20)] = rand(rng, collect("ACGT"))
    end
    svs = [vtest_sv(r, kind, pos, stop, "sv")]
    for j in 1:rand(rng, 0:3)
        q = rand(rng, 30:(length(r) - 30))
        edit = vtest_edit(r, q, q + 1, vtest_alternative(r[q]), "snp$j")
        push!(svs, (edit, vtest_record(edit), vtest_key(edit, "chr1")))
    end
    unique!(s -> s[1].pos, svs)
    return r, svs, [LongDNA{4}(String(guide))], d
end
