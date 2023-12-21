struct BinaryFuseFilterDB
    dbi::DBInfo
    mpt::PathTemplates
    ambig::Union{AmbigIdx, Nothing}
end


function restrictDistance(bffDB::BinaryFuseFilterDB, distance::Int)
    mpt = restrictDistance(bffDB.mpt, distance)
    return BinaryFuseFilterDB(bffDB.dbi, mpt, bffDB.ambig)
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
    precision::DataType = UInt16)
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
    precision::DataType = UInt32)

    dbi = DBInfo(genomepath, name, motif)
    @info "Building Motif templates..."
    mpt = build_PathTemplates(motif)

    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)
    ambig = Vector{LongDNA{4}}()

    for chrom_name in dbi.gi.chrom
        record = reader[chrom_name] # this is possible only with index!
        @info "Working on $chrom_name"
        chrom = dbi.gi.is_fa ? FASTA.sequence(LongDNA{4}, record) : TwoBit.sequence(LongDNA{4}, record)
        guides_fwd = Vector{UInt64}()
        pushguides!(guides_fwd, ambig, dbi, chrom, false; remove_pam = false) # we keep PAM seqeunces - this creates more loci -> more precise filter
        guides_fwd = unique(guides_fwd)
        bff_fwd = BinaryFuseFilter{precision}(guides_fwd; seed = seed, max_iterations = max_iterations)
        if (!all(in.(guides_fwd, Ref(bff_fwd)))) 
            throw("Not all guides are inside the Binary Fuse Filter... Report to the developers.")
        end
        guides_rve = Vector{UInt64}()
        pushguides!(guides_rve, ambig, dbi, chrom, true; remove_pam = false) # guides here will be GGN...EXT and TTN...EXT
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

    save(BinaryFuseFilterDB(dbi, mpt, ambig), 
        joinpath(storage_path, "BinaryFuseFilterDB.bin"))

    @info "Finished."
    return 
end


function not_duplicated_and_in_db(x::Vector{UInt64}, bff::BinaryFuseFilter)
    s = Set(Vector{UInt64}())
    b = BitVector(zeros(length(x)))
    for (i, xi) in enumerate(x)
        if !(xi in s)
            push!(s, xi)
            b[i] = xi in bff
        end
    end
    return b
end


function search_chrom2(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    bffddbir::String, 
    fmidbdir::String, 
    bffDB::BinaryFuseFilterDB)

    guides_uint64 = guide_to_template_format.(copy(guides); alphabet = ALPHABET_TWOBIT)
    guides_uint64_rc = guide_to_template_format.(copy(guides); alphabet = ALPHABET_TWOBIT) 
    guides_fmi = guide_to_template_format.(copy(guides); alphabet = ALPHABET_UINT8)
    guides_fmi_rc = guide_to_template_format.(copy(guides), true; alphabet = ALPHABET_UINT8) # complements guide GGN -> CCN, TTTN -> AAAN

    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = getchromseq(gi.is_fa, reader[chrom])
    fmi = load(joinpath(fmidbdir, chrom * ".bin"))
    bff = load(joinpath(bffddbir, "BinaryFuseFilterDB_" * chrom * ".bin")) # bff_fwd
    detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
    detail_file = open(detail_path, "w")

    # wroking on this guide and his all possible off-targets
    for (i, g) in enumerate(guides) # for each guide
        if bffDB.mpt.motif.extends5
            guide_stranded = reverse(g)
        else
            guide_stranded = g
        end

        # STEP 1. Check in hash whether this OT is there or not
        ot_uint64 = guides_uint64[i][bffDB.mpt.paths] # with PAM PAMseqEXT
        ot_uint64 = map(asUInt64, eachrow(ot_uint64))
        ot_uint64_rc = guides_uint64_rc[i][bffDB.mpt.paths] # with PAMseqEXT - normalized always
        ot_uint64_rc = map(asUInt64, eachrow(ot_uint64_rc))
        # further reduce non-unique seqeunces
        ot_uint64 = not_duplicated_and_in_db(ot_uint64, bff.bff_fwd) # these are both unique and in HashDB
        ot_uint64_rc = not_duplicated_and_in_db(ot_uint64_rc, bff.bff_rve) # BitVector relative to our paths vector

        # STEP 2. actually find the location of the OTs in the genome
        ot = guides_fmi[i][bffDB.mpt.paths[ot_uint64, :]] # GGN + 20N + extension
        ot_rc = guides_fmi_rc[i][bffDB.mpt.paths[ot_uint64_rc, :]] # CCN + 20N + extension
        distances = bffDB.mpt.distances[ot_uint64]
        distances_rc = bffDB.mpt.distances[ot_uint64_rc]
        if bffDB.mpt.motif.extends5
            reverse!(ot; dims = 2) # extension + 20N + NGG
        else # Cpf1
            reverse!(ot_rc; dims = 2) # extension + 20N + NAAA
        end
        
        ot_len = size(ot)[2]
        for i in 1:size(ot)[1] # for each potential OT in fwd
            ot_i = @view ot[i, :]

            fwd_iter = locate(ot_i, fmi)
            if bffDB.mpt.motif.extends5 
                fwd_iter = fwd_iter .+ ot_len .- 1
            end

            for pos in fwd_iter
                line = string(guide_stranded) * "," * "no_alignment" * "," * 
                    string(LongDNA{4}(reinterpret(DNA, ot_i))) * "," * string(distances[i]) * "," *
                    chrom * "," * string(pos) * "," * "+" * "\n"
                write(detail_file, line)
            end
        end

        for i in 1:size(ot_rc)[1] # for each potential OT in rve
            ot_rc_i = @view ot_rc[i, :]
                
            rve_iter = locate(ot_rc_i, fmi)
            if !bffDB.mpt.motif.extends5
                rve_iter = rve_iter .+ ot_len .- 1
            end

            for pos in rve_iter
                line = string(guide_stranded) * "," * "no_alignment" * "," * 
                    string(LongDNA{4}(reinterpret(DNA, ot_rc_i))) * "," * string(distances_rc[i]) * "," *
                    chrom * "," * string(pos) * "," * "-" * "\n"
                write(detail_file, line)
            end
        end
    end
    close(detail_file)
    
    return 
end



"""
```
earch_binaryFuseFilterDB(
    bffddbir::String, 
    fmidbdir::String,
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3)
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
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 0)

    if any(isambig.(guides)) # TODO?
        throw("Ambiguous bases are not allowed in guide queries.")
    end

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
    ThreadsX.map(ch -> search_chrom2(ch, dirname(output_file), guides_, bffddbir, fmidbdir, bffDB), bffDB.dbi.gi.chrom)
    cleanup_detail(output_file)
    return 
end