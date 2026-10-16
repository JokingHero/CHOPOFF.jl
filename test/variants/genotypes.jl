using Random
isdefined(@__MODULE__, :vtest_text) || include(joinpath(@__DIR__, "..", "helpers", "vcf.jl"))
isdefined(@__MODULE__, :vtest_check_tokens) || include(joinpath(@__DIR__, "..", "helpers", "vcf_genotypes.jl"))

@testset "VCF genotypes and haplotype tokens [D01]" begin
    @testset "GT parsing" begin
        gt(s) = (a = UInt8[]; p = CHOPOFF.parse_variant_gt!(a, s); (Int.(a), p))
        m = Int(CHOPOFF.VARIANT_GT_MISSING)
        @test gt("0|1") == ([0, 1], true)
        @test gt("0/1") == ([0, 1], false)
        @test gt("1") == ([1], true)
        @test gt("|0|1") == ([0, 1], true)
        @test gt("/0|1") == ([0, 1], false)
        @test gt("0|1/2") == ([0, 1, 2], false)
        @test gt(".|1") == ([m, 1], true)
        @test gt("./.") == ([m, m], false)
        @test gt("1|2|0") == ([1, 2, 0], true)
        for bad in ("0|", "a/1", "", "0||1", "255/0")
            @test_throws ErrorException gt(bad)
        end
        @test CHOPOFF.variant_subfield("0|1:12:x", 2) == "12"
        @test CHOPOFF.variant_subfield("0|1", 2) === nothing
    end

    # Records: (source, group label, GT string). Labels are :chrom, :unphased or
    # a PS integer; the oracle reads them, the store reads only the GT text.
    function store_from(samples, records)
        store = CHOPOFF.ScanGenotypes(samples)
        buffer = UInt8[]
        for (source, values) in records
            text = ".\tPASS\t.\tGT:PS\t" * join(values, '\t')
            CHOPOFF.push_variant_genotypes!(store, buffer, source, text, eachindex(samples), "test:$source")
        end
        return store
    end
    gt_text(alleles, label) = join((a < 0 ? "." : string(a) for a in alleles),
        label === :unphased ? '/' : '|') * ":" * (label isa Integer ? string(label) : ".")

    @testset "Store layout" begin
        store = store_from(["a", "b%c"], [(3, ["0|0:.", "0/1:."]), (5, ["0|0:.", "0|0:."]), (7, ["1|0:40", ".|0:."])])
        @test store.record_source == [3, 7]
        @test CHOPOFF.variant_record_entries(store, 5) == 1:0
        @test store.sample[CHOPOFF.variant_record_entries(store, 7)] == [1, 2]
        @test store.group == [-3, 41, 0]
        @test CHOPOFF.haplotype_tokens(store, [(7, 1)], []) == ("a:1@40;b%25c:?", 1, 1)
        @test CHOPOFF.haplotype_tokens(store, [], [(3, 1)]; lost = true) == ("b%25c:*", 1, 0)
        @test CHOPOFF.haplotype_tokens(store, [], [(5, 1)]; lost = true) == ("", 0, 0)
        @test_throws ErrorException CHOPOFF.haplotype_tokens(store, [], [(3, 1)])
    end

    # Brute force over assignments: a permutation of every phase group relative
    # to a frame group, and a concrete allele for every missing slot.
    function oracle(records, ploidy, required, forbidden, lost)
        labels = unique(r.label === :unphased || ploidy == 1 ? (:u, r.source) : r.label for r in records)
        label(r) = r.label === :unphased || ploidy == 1 ? (:u, r.source) : r.label
        missing_slots = [(i, j) for (i, r) in enumerate(records) for j in 1:ploidy if r.alleles[j] < 0]
        perms = collect(Iterators.product((collect(permutations_of(ploidy)) for _ in labels)...))
        fills = collect(Iterators.product((0:2 for _ in missing_slots)...))
        # truth[frame][assignment] = copies (frame slot indexes) that hold
        truth = Dict{Any, Vector{Vector{Int}}}()
        for (fi, frame) in enumerate(labels)
            truth[frame] = Vector{Int}[]
            for perm in perms
                perm[fi] == Tuple(1:ploidy) || continue
                for fill in fills
                    allele(i, slot) = (k = findfirst(==((i, slot)), missing_slots);
                        k === nothing ? records[i].alleles[slot] : fill[k])
                    holds = Int[]
                    for j in 1:ploidy
                        ok = !lost
                        for (i, r) in enumerate(records)
                            a = allele(i, perm[findfirst(==(label(r)), labels)][j])
                            if lost
                                ok |= a in forbidden[r.source]
                            else
                                ok &= get(required, r.source, a) == a && !(a in forbidden[r.source])
                            end
                        end
                        ok && push!(holds, j)
                    end
                    push!(truth[frame], holds)
                end
            end
        end
        return labels, truth
    end
    permutations_of(n) = n == 1 ? [(1,)] : n == 2 ? [(1, 2), (2, 1)] :
        [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]

    @testset "Tokens match assignment enumeration" begin
        rng = MersenneTwister(20261014)
        checked = 0
        for case in 1:3000
            ploidy = rand(rng, (1, 2, 2, 2, 3))
            n = rand(rng, 1:4)
            choices = (:chrom, :unphased, 10, 20)
            records = [(; source = 10i, label = rand(rng, choices),
                alleles = [rand(rng, (-1, 0, 0, 1, 1, 2)) for _ in 1:ploidy]) for i in 1:n]
            lost = rand(rng, Bool)
            required = Dict{Int, Int}()
            forbidden = Dict(r.source => Int[] for r in records)
            for r in records
                role = rand(rng, 1:3)
                !lost && role == 1 && (required[r.source] = rand(rng, 1:2))
                (lost || role >= 2) && append!(forbidden[r.source], rand(rng, ([1], [2], [1, 2])))
            end
            !lost && isempty(required) && (required[records[1].source] = 1)
            # A row never both takes and excludes one allele.
            for (s, a) in required
                filter!(!=(a), forbidden[s])
            end
            store = store_from(["s"], [(r.source, [gt_text(r.alleles, r.label)]) for r in records])
            taken = [(s, a) for (s, a) in required]
            excluded = [(s, a) for (s, as) in forbidden for a in as]
            tokens, certain, possible = CHOPOFF.haplotype_tokens(store, taken, excluded; lost)
            # Hom-ref records carry no entry; the oracle sees them as all-0 copies.
            labels, truth = oracle(records, ploidy, required, forbidden, lost)
            any_holds = any(!isempty(h) for t in values(truth) for h in t)
            @test isempty(tokens) == !any_holds
            isempty(tokens) && continue
            parts = [last(rsplit(t, ':'; limit = 2)) for t in split(tokens, ';')]
            ok = vtest_check_tokens(parts, truth, ploidy; lost)
            @test ok
            ok || @info "token mismatch" records=repr(records) ploidy required forbidden lost tokens
            @test certain == count(!=("?"), parts) && possible == ("?" in parts)
            checked += 1
        end
        @test checked > 1000
    end

    @testset "Batch parsing matches sequential parsing" begin
        rng = MersenneTwister(3)
        samples = ["a", "b", "c", "d"]
        values = ["0|0", "0/0", "0|1:5", "1|0:7", "1/1:.", "./.:.", ".|1:5", "1", "0/2:.", "2|1:7"]
        lines = map(1:300) do i
            rand(rng) < 0.05 && return ""
            chrom = rand(rng) < 0.2 ? "chrY" : "chrX"
            "$chrom\t$i\t.\tA\tC,G\t.\t.\t.\tGT:PS\t" * join(rand(rng, values, 4), '\t')
        end
        chrom_index = Dict("chrX" => 1, "chrY" => 2)
        selected = [true, false]
        columns = [1, 3, 4]
        sources = cumsum([!isempty(l) for l in lines])
        sequential = CHOPOFF.ScanGenotypes(samples)
        expected = zeros(UInt128, length(lines))
        buffer = UInt8[]
        for (i, line) in enumerate(lines)
            startswith(line, "chrX") || continue
            fields = split(line, '\t'; limit = 6)
            expected[i] = CHOPOFF.push_variant_genotypes!(sequential, buffer, sources[i], fields[6], columns, "x")
        end
        batched = CHOPOFF.ScanGenotypes(samples)
        carried = vcat((CHOPOFF.push_variant_genotype_batch!(batched, lines[r], sources[r], chrom_index, selected, columns)
            for r in (1:70, 71:71, 72:300))...)
        @test carried == expected
        for field in fieldnames(CHOPOFF.ScanGenotypes)
            @test getfield(batched, field) == getfield(sequential, field)
        end
    end
end
