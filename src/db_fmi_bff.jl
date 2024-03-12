struct BinaryFuseFilterDB
    dbi::DBInfo
    mpt::PathTemplates
    ambig::Union{AmbigIdx, Nothing}
    restrict_to_len::Union{Int, Nothing}
end


function restrictDistance(bffDB::BinaryFuseFilterDB, distance::Int)
    mpt = restrictDistance(bffDB.mpt, distance)
    return BinaryFuseFilterDB(bffDB.dbi, mpt, bffDB.ambig, bffDB.restrict_to_len)
end


struct BinaryFuseFilterDBperChrom{K<:Union{UInt8, UInt16, UInt32}}
    dbi::DBInfo
    bff_fwd::BinaryFuseFilter{K}
    bff_rve::BinaryFuseFilter{K}
    chrom::String
end


"""
```
build_binaryFuseFilterDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storage_dir::String;
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    precision::DataType = UInt32)
```

Prepare hashDB index for future searches using `search_hashDB`.


# Arguments
`name` - Your preferred name for this index to ease future identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for and at what maxium distance.

`storage_dir`  - Directory to the where many files needed by the database will be saved. Naming 
                 of the files follows this pattern: 
                 BinaryFuseFilterDB_ + chromsome + .bin
                 Each unique motif has its own file naming created.

`seed`  - Optional. Seed is used during hashing for randomization.

`max_iterations` - When finding hashing structure for binary fuse filter it might fail sometimes, 
                   we will retry `max_iterations` number of times though.
                
`precision`- The higher the precision the larger the database, but also chances for error decrease dramatically.
             We support UInt8, UInt16, and UInt32.

`restrict_to_len` - Restrict lengths of the `motif` for the purpose of checking its presence in the genome.
                    Allows for significant speedups when expanding all possible sequences for each guide, as we will expand
                    up to the specified length here. For example, default setting for Cas9, would restrict standard 20bp to
                    16bp for the genome presence check, for distance of 4 that would be 8 bases (4bp from the 20 - 16, and 4 
                    bases because of the potential extension) that have to be actually aligned in the genome.

# Examples
```julia-repl
$(make_example_doc("binaryFuseFilterDB"))
```
"""
function build_binaryFuseFilterDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storage_path::String;
    seed::UInt64 = UInt64(0x726b2b9d438b9d4d),
    max_iterations::Int = 10,
    precision::DataType = UInt32,
    restrict_to_len::Int = length_noPAM(motif) - motif.distance)

    if restrict_to_len > (length_noPAM(motif) + motif.distance)
        restrict_to_len = length_noPAM(motif) + motif.distance
        @warn "Removing length restriction, expect this to be slow and possibly explode your memory!"
    end

    dbi = DBInfo(genomepath, name, motif)
    @info "Building Motif templates..."
    # PAM here adds too many sequences to verify...
    mpt = build_PathTemplates(motif; restrict_to_len = restrict_to_len, withPAM = false)

    mkpath(storage_path)
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    ambig = Vector{LongDNA{4}}() # TODO

    for chrom_name in dbi.gi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.gi.is_fa ? FASTA.sequence(LongDNA{4}, record) : TwoBit.sequence(LongDNA{4}, record)
        guides_fwd = Vector{UInt64}()
        pushguides!(guides_fwd, ambig, dbi, chrom, false; remove_pam = true, restrict_to_len = restrict_to_len) # we need to check PAM on the genome
        guides_fwd = unique(guides_fwd)
        bff_fwd = BinaryFuseFilter{precision}(guides_fwd; seed = seed, max_iterations = max_iterations)
        if (!all(in.(guides_fwd, Ref(bff_fwd)))) 
            throw("Not all guides are inside the Binary Fuse Filter... Report to the developers.")
        end
        guides_rve = Vector{UInt64}()
        pushguides!(guides_rve, ambig, dbi, chrom, true; remove_pam = true, restrict_to_len = restrict_to_len) # guides here will be GGN...EXT and TTN...EXT
        guides_rve = unique(guides_rve)
        bff_rve = BinaryFuseFilter{precision}(guides_rve; seed = seed, max_iterations = max_iterations)
        if (!all(in.(guides_rve, Ref(bff_rve)))) 
            throw("Not all guides are inside the Binary Fuse Filter... Report to the developers.")
        end
        save(BinaryFuseFilterDBperChrom(dbi, bff_fwd, bff_rve, chrom_name), 
            joinpath(storage_path, "BinaryFuseFilterDB_" * chrom_name * ".bin"))
    end

    # TODO add ambiguity handling here?!
    ambig = nothing # ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing
    close(ref)

    save(BinaryFuseFilterDB(dbi, mpt, ambig, restrict_to_len), 
        joinpath(storage_path, "BinaryFuseFilterDB.bin"))

    @info "Finished."
    return 
end


function search_chrom2(
    chrom_name::String,
    detail::String, 
    guides::Vector{LongDNA{4}},
    bffddbir::String, 
    fmidbdir::String, 
    genomepath::String,
    bffDB::BinaryFuseFilterDB)

    ref = open(genomepath, "r")
    reader = bffDB.dbi.gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    chrom = getchromseq(bffDB.dbi.gi.is_fa, reader[chrom_name])

    guides_uint64 = guide_to_template_format.(copy(guides); alphabet = ALPHABET_TWOBIT)
    guides_uint64_rc = guide_to_template_format.(copy(guides); alphabet = ALPHABET_TWOBIT) 
    guides_fmi = guide_to_template_format.(copy(guides); alphabet = ALPHABET_UINT8)
    guides_fmi_rc = guide_to_template_format.(copy(guides), true; alphabet = ALPHABET_UINT8) # complements guide GGN -> CCN, TTTN -> AAAN
    guides_ambig = guide_to_template_format_ambig.(copy(guides))
    guides_ambig_rc = guide_to_template_format_ambig.(copy(guides), true)

    fmi = load(joinpath(fmidbdir, chrom_name * ".bin"))
    bff = load(joinpath(bffddbir, "BinaryFuseFilterDB_" * chrom_name * ".bin"))
    restrict_to_len = bffDB.mpt.restrict_to_len
    detail_path = joinpath(detail, "detail_" * chrom_name * ".csv")
    detail_file = open(detail_path, "w")

    # wroking on this guide and his all possible off-targets
    pam = bffDB.mpt.motif.fwd[bffDB.mpt.motif.pam_loci_fwd]
    pam_rc = bffDB.mpt.motif.rve[bffDB.mpt.motif.pam_loci_rve]
    pam_len = length(pam)
    tail_len = size(bffDB.mpt.paths)[2] - restrict_to_len # no PAM assumptions
    for (i, g) in enumerate(guides)
        if bffDB.mpt.motif.extends5
            guide_stranded = reverse(g)
        else
            guide_stranded = g
        end

        # STEP 1. Check in hash whether this OT is there or not
        ot_uint64 = guides_uint64[i][bffDB.mpt.paths[:, 1:restrict_to_len]] # PAMseqEXT
        ot_uint64 = map(x -> ARTEMIS.asUInt(UInt64, x), eachrow(ot_uint64))
        ot_uint64_rc = guides_uint64_rc[i][bffDB.mpt.paths[:, 1:restrict_to_len]] # PAMseqEXT - normalized always
        ot_uint64_rc = map(x -> asUInt(UInt64, x), eachrow(ot_uint64_rc))
        # further reduce non-unique seqeunces
        ot_uint64 = in.(ot_uint64, Ref(bff.bff_fwd))
        ot_uint64_rc = in.(ot_uint64_rc, Ref(bff.bff_rve))

        # STEP 2. actually find the location of the OTs in the genome
        ot = guides_fmi[i][bffDB.mpt.paths[ot_uint64, 1:restrict_to_len]] # GGN + 20N + extension
        ot_rc = guides_fmi_rc[i][bffDB.mpt.paths[ot_uint64_rc, 1:restrict_to_len]] # CCN + 20N + extension
        ot_tail = guides_ambig[i][bffDB.mpt.paths[ot_uint64, restrict_to_len + 1:end]] # GGN + 20N + extension
        ot_rc_tail = guides_ambig_rc[i][bffDB.mpt.paths[ot_uint64_rc, restrict_to_len + 1:end]] # CCN + 20N + extension

        distances = bffDB.mpt.distances[ot_uint64]
        distances_rc = bffDB.mpt.distances[ot_uint64_rc]
        if bffDB.mpt.motif.extends5
            reverse!(ot; dims = 2) # extension + 20N + NGG
            reverse!(ot_tail; dims = 2)
        else # Cas12a
            reverse!(ot_rc; dims = 2) # extension + 20N + NAAA
            reverse!(ot_rc_tail; dims = 2)
        end
        

        for j in 1:size(ot)[1]
            ot_j = @view ot[j, :]

            fwd_iter = ARTEMIS.locate(ot_j, fmi)
            for pos in fwd_iter

                if bffDB.mpt.motif.extends5
                    pass = all(ARTEMIS.iscompatible.(
                        pam, chrom[(pos + restrict_to_len):(pos + restrict_to_len + pam_len - 1)])) && 
                        all(ARTEMIS.iscompatible.(ot_tail[j, :], chrom[(pos - tail_len):(pos - 1)]))
                else
                    pass = all(ARTEMIS.iscompatible.(
                        pam, chrom[(pos - pam_len):(pos - 1)])) && all(ARTEMIS.iscompatible.(
                        ot_tail[j, :], chrom[(pos + restrict_to_len):(pos + restrict_to_len + tail_len - 1)]))
                end

                if pass
                    # TODO remove trailing Ns?!
                    if bffDB.mpt.motif.extends5
                        aln_ref = chrom[(pos - tail_len):pos + restrict_to_len - 1]
                        pos = pos + restrict_to_len + tail_len - 1 - 1
                    else
                        aln_ref = chrom[pos:pos + restrict_to_len + tail_len - 1]
                        pos = pos - pam_len
                    end

                    line = string(guide_stranded) * "," * "no_alignment" * "," * 
                        string(aln_ref) * "," * string(distances[j]) * "," *
                        chrom_name * "," * string(pos) * "," * "+" * "\n"
                    write(detail_file, line)
                end
            end
        end

        for j in 1:size(ot_rc)[1]
            ot_rc_j = @view ot_rc[j, :]
                
            rve_iter = ARTEMIS.locate(ot_rc_j, fmi)
            for pos in rve_iter

                if bffDB.mpt.motif.extends5
                    pass = all(ARTEMIS.iscompatible.(
                        pam_rc, chrom[(pos - pam_len):(pos - 1)])) && all(ARTEMIS.iscompatible.(
                        ot_rc_tail[j, :], chrom[(pos + restrict_to_len):(pos + restrict_to_len + tail_len - 1)]))
                else
                    pass = all(ARTEMIS.iscompatible.(
                        pam_rc, chrom[(pos + restrict_to_len):(pos + restrict_to_len + pam_len - 1)])) && 
                        all(ARTEMIS.iscompatible.(ot_rc_tail[j, :], chrom[(pos - tail_len):(pos - 1)]))
                end

                if pass
                    if bffDB.mpt.motif.extends5
                        aln_ref = chrom[pos:pos + restrict_to_len + tail_len - 1]
                        pos = pos - pam_len
                    else
                        aln_ref = chrom[(pos - tail_len):pos + restrict_to_len - 1]
                        pos = pos + restrict_to_len + tail_len - 1
                    end
                    
                    line = string(guide_stranded) * "," * "no_alignment" * "," * 
                        string(aln_ref) * "," * string(distances_rc[j]) * "," *
                        chrom_name * "," * string(pos) * "," * "-" * "\n"
                    write(detail_file, line)
                end
            end
        end
    end

    close(ref)
    close(detail_file)
    
    return 
end



"""
```
search_binaryFuseFilterDB(
    bffddbir::String, 
    fmidbdir::String,
    genomepath::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 0)
```

Find all off-targets for `guides` within distance of `dist` using BinaryFuseFilterDB located at `storage_dir`.

Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

```
5' - ...ACGTCATCG NGG - 3'  -> will be input: ...ACGTCATCG
     guide        PAM
    
3' - CCN GGGCATGCT... - 5'  -> will be input: ...AGCATGCCC
     PAM guide
```

# Arguments

`bffdbdir` - Folder location where BinaryFuseFilterDB is stored at.

`fmidbdir` - Folder location where FM-index is build.

`guides` - Vector of your guides, without PAM.

`output_file` - Path and name for the output file, this will be comma separated table, therefore `.csv` extension is preferred. 
This search will create intermediate files which will have same name as `output_file`, but with a sequence prefix. Final file
will contain all those intermediate files.

`distance` - Defines maximum levenshtein distance (insertions, deletions, mismatches) for 
which off-targets are considered.

# Examples
```julia-repl
$(make_example_doc("binaryFuseFilterDB"))
```
"""
function search_binaryFuseFilterDB(
    bffddbir::String, 
    fmidbdir::String,
    genomepath::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 0)

    if any(isambig.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    #gi2 = load(joinpath(fmidbdir, "genomeInfo.bin"))
    #gi = GenomeInfo(genomepath)
    #if !isequal(gi, gi2)
    #    msg = "Supplied genome is different than the genome used for building of FM-index. "
    #    if isequal(Set(gi.chrom), Set(gi.chrom))
    #        msg *= "Supplied genome has the same chromosome names."
    #    else
    #        if all(occursin.(gi.chrom, gi2.chrom))
    #            msg *= "Supplied genome has different chromosome names, but all exist in FM-index."
    #        else
    #            throw("Supplied genome has different chromosome names, but not all exist in FM-index.")
    #        end
    #    end
    #    @warn msg
    #end

    guides_ = copy(guides)
    bffDB = load(joinpath(bffddbir, "BinaryFuseFilterDB.bin"))

    if any(length_noPAM(bffDB.dbi.motif) .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end

    # reverse guides so that PAM is always on the left
    if bffDB.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    bffDB = restrictDistance(bffDB, distance)
    # we input guides that are in forward search configuration # which means 20bp-NGG
    ThreadsX.map(ch -> search_chrom2(ch, dirname(output_file), guides_, bffddbir, fmidbdir, genomepath, bffDB), bffDB.dbi.gi.chrom)
    
    cleanup_detail(output_file)
    return 
end