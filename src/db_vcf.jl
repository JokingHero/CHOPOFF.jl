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


"""
```
build_vcfDB(
    name::String, 
    genomepath::String, 
    vcfpath::String,
    motif::Motif,
    storage_path::String,
    hash_len::Int = min(length_noPAM(motif) - motif.distance, 16))
```

Builds a database of all potential off-targets that overlap any of the variants in the VCF file.
It supports combinations of variants that are close to each other, will report all possible combinations of 
variants. This database uses simialr principles to `prefixHashDB`, also utilizes hashed prefix of specific length.


# Arguments
`name` - Your preferred name for this index for future identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`vcfpath` - Path to the VCF file, it has to be compatible with your genome.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_path`  - Path to the where index will be saved. Normally, it includes ".bin" extension.

`hash_len` - length of the prefix that is stored inside the hash


# Examples
```julia
$(make_vcf_example_doc())
```
"""
function build_vcfDB(
    name::String, 
    genomepath::String, 
    vcfpath::String,
    motif::Motif,
    storage_path::String,
    hash_len::Int = min(length_noPAM(motif) - motif.distance, 16))

    dbi = DBInfo(genomepath, name, motif; vcf_filepath = vcfpath)
    hash_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    suffix_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", (CHOPOFF.length_noPAM(motif) - hash_len + motif.distance) * 2); base = 2))
    ot_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", (CHOPOFF.length_noPAM(motif) + motif.distance) * 2); base = 2))

    rs_ids, rs_chroms, rs_ref, rs_ranges, rs_alt = CHOPOFF.parse_vcf(vcfpath)
    l = length_noPAM(motif) + motif.distance
    lp = length(motif) + motif.distance # with PAM

    @info "Step 1: Parsing the genomic relation to the VCF file."
        # For each chromosome parallelized we build database
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    ambig_guides = Vector{LongDNA{4}}()
    ambig_chrom = Vector{dbi.gi.chrom_type}()
    ambig_pos = Vector{dbi.gi.pos_type}()
    ambig_isplus = Vector{Bool}()
    ambig_annot = Vector{String}() # change to InlineStrings

    for ch in unique(rs_chroms)
        #ch = first(unique(rs_chroms)) # REMOVE
        ch_ = convert(dbi.gi.chrom_type, findfirst(isequal(ch), dbi.gi.chrom))

        chrom_seq = CHOPOFF.getchromseq(dbi.gi.is_fa, reader[ch])
        chrom_idx = findall(isequal(ch), rs_chroms)
        # now for every snp (we assume they are sorted)
        # group snps by proximity - as we have to enumerate all permutations of rs_alt for them
        grouping = 1
        grouping_idx = ones(Int, length(chrom_idx))
        for i in 1:(length(chrom_idx) - 1)
            x = (rs_ranges[chrom_idx[i]].start - l):(rs_ranges[chrom_idx[i]].stop + l)
            if length(intersect(x, rs_ranges[chrom_idx[i+1]])) > 0 # rs overlaps
                grouping_idx[i + 1] = grouping
            else
                grouping += 1
                grouping_idx[i + 1] = grouping
            end
        end

        # for each group we analyze these snps together
        for group in unique(grouping_idx)
            # group = first(unique(grouping_idx))
            idxs = chrom_idx[grouping_idx .== group]
            if length(idxs) == 1 # simple case - singular snp - potentially many alternate alleles
                idxs = idxs[1]
                # -----A-----
                #      5     # SNP start
                # --2-----8--# lp of 4bp 
                seq_start = rs_ranges[idxs].start - lp + 1
                seq_end = rs_ranges[idxs].stop + lp - 1
                seq = chrom_seq[seq_start:seq_end] # lp * 2 - 1 (one base is shared between lp)
                # annot in style of rsID:REF/ALT
                seq_ref = string(seq[lp:lp + length(rs_ranges[idxs]) - 1])
                if seq_ref != rs_ref[idxs]
                    @warn "Possible mismatch between VCF REF column:" * 
                        rs_ids[idxs] * ";" * ch * ":" * string(rs_ranges[idxs].start) *
                        string(rs_ranges[idxs].stop) * " and the reference file."
                end
                for alt in rs_alt[idxs]
                    seq_alt = copy(seq)
                    seq_alt = seq_alt[1:(lp - 1)] * alt * seq_alt[lp + length(seq_ref):end]
                    g, g_rev, gp, gp_rev = CHOPOFF.gatherofftargets(seq_alt, dbi)
                    overlaps = in.(lp, gp) .| in.(lp + length(alt) - 1, gp)
                    g = g[overlaps]
                    gp = gp[overlaps]
                    gp = Base.map(x -> seq_start + x.stop, gp)
                    overlaps = in.(lp, gp_rev) .| in.(lp + length(alt) - 1, gp_rev)
                    g_rev = g_rev[overlaps]
                    gp_rev = gp_rev[overlaps]
                    gp_rev = Base.map(x -> seq_start + x.start - 1, gp_rev)
                    this_annot = join(rs_ids[idxs], ";") * ":" *seq_ref * "/" * string(alt)
                    append!(ambig_guides, g)
                    append!(ambig_guides, g_rev)
                    append!(ambig_chrom, repeat([ch_], length(g) + length(g_rev)))
                    append!(ambig_pos, gp)
                    append!(ambig_pos, gp_rev)
                    append!(ambig_isplus, trues(length(g)))
                    append!(ambig_isplus, falses(length(g_rev)))
                    append!(ambig_annot, repeat([this_annot], length(g) + length(g_rev)))
                end
            else # complex case multiple overlapping snps, can have multiple alternate alleles
                alt_combs = Iterators.product(rs_alt[idxs]...)
                # we need to track which regions belong to which snp
                # and which guide overlaps which snps - 
                # we allow unique guides by sequence and annotations
                for alt_comb in alt_combs
                    seq_alt = dna""
                    rg_comb = Vector{UnitRange{Int64}}()
                    annot_comb = Vector{String}()
                    seq_start = rs_ranges[idxs][1].start - lp + 1
                    for (i, alt) in enumerate(alt_comb)
                        rst = rs_ranges[idxs[i]].start
                        rsp = rs_ranges[idxs[i]].stop
                        alt_len = length(alt) == 0 ? -1 : length(alt) # deletions 
                        #----ref--ref-------
                        #  --alt--alt--
                        #  --rng--rng-- range of the alt on the alt comb
                        seq_ref = string(chrom_seq[rst:rsp])
                        if seq_ref != rs_ref[idxs[i]]
                            @warn "Possible mismatch between VCF REF column:" * 
                            rs_ids[idxs[i]] * ";" * ch * ":" * string(rs_ranges[idxs[i]].start) *
                            string(rs_ranges[idxs[i]].stop) * " and the reference file."
                        end
                        push!(annot_comb, join(rs_ids[idxs][i], ";") * ":" *seq_ref * "/" * string(alt))
                        if i == 1 # first snp --alt
                            seq_alt *= chrom_seq[(rst - lp + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                        elseif i == length(alt_comb) # last snp --alt--
                            # first --alt part
                            rsp_prev = rs_ranges[idxs[i-1]].stop
                            seq_alt *= chrom_seq[(rsp_prev + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                            # and -- ending part
                            seq_alt *= chrom_seq[(rsp + 1):(rsp + lp - 1)]
                        else # midle snps --alt part
                            rsp_prev = rs_ranges[idxs[i-1]].stop
                            seq_alt *= chrom_seq[(rsp_prev + 1):(rst - 1)] * alt
                            push!(rg_comb, (length(seq_alt) - alt_len + 1):length(seq_alt))
                        end
                    end
                    # guides here are without PAM, but with extensions, all in GGN-EXT config
                    # guide pos here have original guide pattern locations (with PAM)
                    g, g_rev, gp, gp_rev = CHOPOFF.gatherofftargets(seq_alt, dbi)
                    for (i, gr) in enumerate(gp)
                        # overlaps for each of the SNPs
                        overlaps = BitVector(Base.map(x -> length(intersect(gr, x)) > 0, rg_comb))
                        if sum(overlaps) > 0
                            push!(ambig_guides, g[i])
                            push!(ambig_chrom, ch_)
                            push!(ambig_pos, seq_start + gp[i].stop)
                            push!(ambig_isplus, true)
                            push!(ambig_annot, join(annot_comb[overlaps], "-"))
                        end
                    end
                    for (i, gr) in enumerate(gp_rev)
                        overlaps = BitVector(Base.map(x -> length(intersect(gr, x)) > 0, rg_comb))
                        if sum(overlaps) > 0
                            push!(ambig_guides, g_rev[i])
                            push!(ambig_chrom, ch_)
                            push!(ambig_pos, seq_start + gp_rev[i].start - 1)
                            push!(ambig_isplus, false)
                            push!(ambig_annot, join(annot_comb[overlaps], "-"))
                        end
                    end
                end
            end
        end
    end
    close(ref)

    @info "Step 2: Constructing Paths for hashes"
    mpt = build_PathTemplates(motif; restrict_to_len = hash_len, withPAM = false)
    paths = mpt.paths[:, 1:hash_len]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto)))) # how can there be no duplicated function?!

    if length(ambig_guides) > 0
        @info "Step 3: Constructing DB for ambigous gRNAs."
        order = sortperm(ambig_guides)
        ambig_guides = ambig_guides[order]
        ambig_chrom = ambig_chrom[order]
        ambig_pos = ambig_pos[order]
        ambig_isplus = ambig_isplus[order]
        ambig_annot = InlineStrings.inlinestrings(ambig_annot[order])    
        save(CHOPOFF.build_ambigPrefixHashDB(
                ambig_guides, ambig_chrom, ambig_pos, ambig_isplus,
                l, hash_len, ot_type, hash_type, suffix_type, ambig_annot, 
                CHOPOFF.SymbolicAlignments(dbi, paths[not_dups, :], mpt.distances[not_dups], hash_len)), 
            storage_path)
        @info "Finished constructing vcfDB in " * storage_path
    else
        @info "No guides were found. vcfDB is not saved."
    end
    return 
end


"""
```
search_vcfDB(
    storage_path::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    early_stopping::Vector{Int} = Int.(floor.(exp.(0:distance))))
```

Find all off-targets for `guides` within distance of `dist` using vcfDB located at `storage_dir`.
Uses early stopping to stop searching when a guide passes a limit on number of off-targets. This method does not 
keep track of the off-target locations and does not filter overlapping off-targets, therefore it might hit the 
early stopping condition a little earlier than intended. Especially, when variants have multiple ALT and 
multiple variants are overlapping off-targets, this function will report each combination of the overlapping variants.
Each of these combinations will also count towards early stopping condition.

Assumes that your guides do not contain PAM, and are all in the same direction as 
you would order for the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM
    
3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`output_file` - Path and name for the output file, this will be comma separated table, therefore `.csv` extension is preferred. 
This search will create intermediate files which will have same name as `output_file`, but with a sequence prefix. Final file
will contain all those intermediate files.

`distance` - Defines maximum levenshtein distance (insertions, deletions, mismatches) for 
which off-targets are considered.

`early_stopping` - Integer vector. Early stopping condition. For example for distance 2, we need vector with 3 values e.g. [1, 1, 5].
Which means we will search with "up to 1 offtargets within distance 0", "up to 1 offtargets within distance 1"...

# Examples
```julia
$(make_vcf_example_doc())
```
"""
function search_vcfDB(
    storage_path::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    early_stopping::Vector{Int} = Int.(floor.(exp.(0:distance))))

    if length(early_stopping) != (distance + 1)
        error("Specify one early stopping condition for a each distance, starting from distance 0.")
    end

    adb = load(storage_path)
    ot_len = CHOPOFF.length_noPAM(adb.mpt.dbi.motif) + adb.mpt.dbi.motif.distance
    ot_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", ot_len * 2); base = 2))
    s_len = ot_len - adb.mpt.hash_len

    if distance > adb.mpt.dbi.motif.distance
        error("For this database maximum distance is " * string(adb.mpt.dbi.motif.distance))
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if adb.mpt.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    paths = adb.mpt.paths[adb.mpt.paths_distances .<= distance, :]
    mkpath(dirname(output_file))

    Base.map(guides_) do g # maybe a function would be faster than lambda here?
        guides_formated = CHOPOFF.guide_to_template_format(g; alphabet = CHOPOFF.ALPHABET_TWOBIT)
        asa = guides_formated[paths]
        asa = Base.map(x -> CHOPOFF.asUInt(eltype(adb.prefix), x), eachrow(asa))
        asa = unique(asa)
        asa = CHOPOFF.potential_ots_idx(asa, adb.prefix)
        
        es_acc = zeros(Int64, length(early_stopping))   
        detail_path = joinpath(dirname(output_file), "detail_" * string(g) * ".csv")
        detail_file = open(detail_path, "w")
        guide_stranded = adb.mpt.dbi.motif.extends5 ? reverse(g) : g
        guide_stranded = string(guide_stranded)
        if length(asa) == 0
            close(detail_file)
            return
        end

        asa = vcat(collect.(asa)...)
        prefixes = adb.prefix[asa]
         asa = adb.prefix_idx[asa] # actual indxes of suffixes
            
        dups = CHOPOFF.duplicated(asa) .& CHOPOFF.duplicated(prefixes)
        asa = asa[.!dups]
        prefixes = prefixes[.!dups]
        suffixes = adb.suffix[asa]
        ots = LongDNA{4}.((convert.(ot_type, prefixes) .<< (2 * s_len)) .| 
                convert.(ot_type, suffixes), ot_len)

        @inbounds for (i, ot) in enumerate(ots)
            idx = asa[i] # actual index for chrom/pos/isplus/annot
            bv_start = (idx - 1) * ot_len + 1
            bv_end = idx * ot_len
            bv = adb.is_ambig[bv_start:bv_end]
            bv_start = sum(adb.is_ambig[1:bv_start]) + 1
            bv_end = bv_start + sum(bv) - 1
            ot[bv] = reinterpret.(DNA, adb.ambig[bv_start:bv_end])
                
            aln = CHOPOFF.align(g, ot, distance, iscompatible)
            if aln.dist <= distance
                if adb.mpt.dbi.motif.extends5
                    aln_guide = reverse(aln.guide)
                    aln_ref = reverse(aln.ref)
                else
                    aln_guide = aln.guide
                    aln_ref = aln.ref
                end
                strand = adb.isplus[idx] ? "+" : "-"
                ot = guide_stranded * "," * aln_guide * "," * 
                    aln_ref * "," * string(aln.dist) * "," *
                    adb.mpt.dbi.gi.chrom[adb.chrom[idx]] * "," * 
                        string(adb.pos[idx]) * "," * strand * "," *
                        string(adb.annot[idx]) * "\n"
                write(detail_file, ot)
                es_acc[aln.dist + 1] += 1
                if es_acc[aln.dist + 1] >= early_stopping[aln.dist + 1]
                    close(detail_file)
                    return
                end
            end
        end
        close(detail_file)
        return
    end

    cleanup_detail(output_file; 
        first_line = "guide,alignment_guide,alignment_reference,distance,chromosome,start,strand,variants\n")
    return
end