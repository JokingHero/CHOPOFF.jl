struct VcfDB
    dbi::DBInfo
    ambig::AmbigIdx
end


function parse_vcf(vcf_filepath::String)
    # https://samtools.github.io/hts-specs/VCFv4.2.pdf
    reader = extension(vcf_filepath) == ".gz" ? 
        VCF.Reader(GzipDecompressorStream(open(vcf_filepath))) : 
        VCF.Reader(open(vcf_filepath))

    rs_ids = Vector{Vector{String}}()
    rs_chroms = Vector{String}()
    rs_ref = Vector{String}()
    rs_ranges = Vector{UnitRange{Int64}}()
    rs_alt = Vector{Vector{LongDNA{4}}}()
    recordNum = 1
    for record in reader
        push!(rs_ids, isempty(record.id) ? Vector([string(recordNum)]) : VCF.id(record))
        push!(rs_chroms, VCF.chrom(record))
        push!(rs_ref, VCF.ref(record))
        push!(rs_ranges, VCF.pos(record):(VCF.pos(record) + length(VCF.ref(record)) - 1))
        if  isempty(record.alt)
            push!(rs_alt, Vector([dna""]))
        else
            r_alt = LongDNA{4}.(VCF.alt(record))
            # compress using ambiguity codes
            r_len = length.(r_alt)
            r_alt_compr = Vector{LongDNA{4}}()
            for l in unique(r_len)
                this_l = r_len .== l
                li = LongDNA{4}()
                for i in 1:l
                    push!(li, TO_AMBIGUOUS[join(sort(getindex.(r_alt[this_l], i)))])
                end
                push!(r_alt_compr, li)
            end
            push!(rs_alt, r_alt_compr)
        end
        recordNum += 1
    end
    close(reader)
    return rs_ids, rs_chroms, rs_ref, rs_ranges, rs_alt
end


function build_vcfDB(
    name::String, 
    genomepath::String, 
    vcfpath::String,
    motif::Motif,
    storagedir::String)

    motif = setambig(setdist(motif, 1), length(motif))
    dbi = DBInfo(genomepath, name, motif; vcf_filepath = vcfpath)
    rs_ids, rs_chroms, rs_ref, rs_ranges, rs_alt = parse_vcf(vcfpath)
    motif_len = length(motif) + motif.distance # include distance in all calculations!

    # For each chromsome paralelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    all_guides = Vector{LongDNA{4}}()
    guide_annot = Vector{String}()

    for ch in unique(rs_chroms)
        #ch = first(unique(rs_chroms)) # REMOVE
        chrom_seq = getchromseq(dbi.gi.is_fa, reader[ch])
        chrom_idx = findall(isequal(ch), rs_chroms)
        # now for every snp (we assume they are sorted)
        # group snps by proximity - as we have to enumerate all permutations of rs_alt for them
        grouping = 1
        grouping_idx = ones(Int, length(chrom_idx))
        for i in 1:(length(chrom_idx) - 1)
            x = (rs_ranges[chrom_idx[i]].start - motif_len):(rs_ranges[chrom_idx[i]].stop + motif_len)
            if length(intersect(x, rs_ranges[chrom_idx[i+1]])) > 0 # rs overlaps
                grouping_idx[i + 1] = grouping
            else
                grouping += 1
                grouping_idx[i + 1] = grouping
            end
        end

        # for each group we analyze these snps together
        for group in unique(grouping_idx)
            idxs = chrom_idx[grouping_idx .== group]
            if length(idxs) == 1 # simple case - singular snp - potentially many alternate alleles
                idxs = idxs[1]
                seq = chrom_seq[(rs_ranges[idxs].start - motif_len + 1):(rs_ranges[idxs].stop + motif_len)]
                guides = Vector{LongDNA{4}}()
                for alt in rs_alt[idxs]
                    seq_alt = copy(seq)
                    seq_alt = seq_alt[1:(motif_len - 1)] * alt * seq_alt[(motif_len + length(rs_ref[idxs])):end - 1]
                    g, gp = gatherofftargets(seq_alt, dbi)
                    append!(guides, g)
                end
                guides = unique(guides)
                append!(all_guides, guides)
                for i in guides
                    push!(guide_annot, join(rs_ids[idxs], ";"))
                end
            else # complex case multiple overlapping snps, can have multiple alternate alleles
                alt_combs = Iterators.product(rs_alt[idxs]...)
                temp_guides = Vector{LongDNA{4}}() # before aplying unique filter!
                temp_annot = Vector{String}()
                # we need to track which regions belong to which snp
                # and which guide overlaps which snps - 
                # we allow unique guides by sequence and annotations
                for alt_comb in alt_combs
                    # alt_comb = first(alt_combs) REMOve me
                    # construct sequence
                    seq_alt = dna""
                    rg_comb = Vector{UnitRange{Int64}}()
                    for (i, alt) in enumerate(alt_comb)
                        rst = rs_ranges[idxs[i]].start
                        rsp = rs_ranges[idxs[i]].stop
                        alt_len = length(alt) == 0 ? -1 : length(alt) # deletions 
                        #----ref--ref-------
                        #  --alt--alt--
                        #  --rng--rng-- range of the alt on the alt comb
                        if i == 1 # first snp --alt
                            seq_alt *= chrom_seq[(rst - motif_len + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                        elseif i == length(alt_comb) # last snp --alt--
                            # first --alt part
                            rsp_prev = rs_ranges[idxs[i-1]].stop
                            seq_alt *= chrom_seq[(rsp_prev + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                            # and -- ending part
                            seq_alt *= chrom_seq[(rsp + 1):(rsp + motif_len)]
                        else # midle snps --alt part
                            rsp_prev = rs_ranges[idxs[i-1]].stop
                            seq_alt *= chrom_seq[(rsp_prev + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                        end
                    end
                    # guides here are without PAM, but with extensions, all in GGN-EXT config
                    # guide pos here have original guide pattern locations (with PAM)
                    g, guide_ranges = gatherofftargets(seq_alt, dbi)
                    for (i, gr) in enumerate(guide_ranges)
                        overlaps = Base.map(x -> length(intersect(gr, x)) > 0, rg_comb)
                        if sum(overlaps) > 0
                            push!(temp_annot, join(join.(rs_ids[idxs[overlaps]], ";"), ";"))
                            push!(temp_guides, g[i])
                        end
                    end
                end
                # now reduce guides by their unique sequence and set of overlaping snps
                temp_merged = map(x -> string(x[1]) * x[2], zip(temp_guides, temp_annot))
                temp_not_dup = map(x -> sum(x .== temp_merged) == 1, temp_merged)
                append!(all_guides, temp_guides[temp_not_dup])
                append!(guide_annot, temp_annot[temp_not_dup])
            end
        end
    end
    close(ref)

    db = VcfDB(dbi, AmbigIdx(all_guides, guide_annot))
    save(db, joinpath(storagedir, "vcfDB.bin"))
    @info "Finished constructing vcfDB in " * storagedir
    @info "Database size is:" *
        "\n length -> " * string(length(db.ambig)) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "vcfDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


function search_vcfDB(
    storagedir::String,
    guides::Vector{LongDNA{4}})

    db = load(joinpath(storagedir, "vcfDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if db.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), 2)
    len_noPAM = length_noPAM(db.dbi.motif)
    for (i, s) in enumerate(guides_)
        res[i, 1] += sum(findbits(s, db.ambig))
        
        d1_combs = LongDNA{4}.(comb_of_d1_extended(string(s))) # 1 distance
        bits_mapped = map(x -> findbits(x, db.ambig), d1_combs)
        res[i, 2] += sum(reduce(.|, bits_mapped))
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, vcat(col_d, :guide))
    return res
end