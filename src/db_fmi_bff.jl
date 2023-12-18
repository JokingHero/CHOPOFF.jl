struct BinaryFuseFilterDB
    dbi::DBInfo
    mpt::PathTemplates
    ambig::Union{AmbigIdx, Nothing}
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
    precision::DataType = UInt16)

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
        pushguides!(guides_fwd, ambig, dbi, chrom, false)
        guides_fwd = unique(guides_fwd)
        bff_fwd = BinaryFuseFilter{precision}(guides_fwd; seed = seed, max_iterations = max_iterations)
        if (!all(in.(guides_fwd, Ref(bff_fwd)))) 
            throw("Not all guides are inside the Binary Fuse Filter... Report to the developers.")
        end
        guides_rve = Vector{UInt64}()
        pushguides!(guides_rve, ambig, dbi, chrom, true)
        guides_rve = unique(guides_rve)
        bff_rve = BinaryFuseFilter{precision}(guides_rve; seed = seed, max_iterations = max_iterations)
        if (!all(in.(guides_rve, Ref(bff_fwd)))) 
            throw("Not all guides are inside the Binary Fuse Filter... Report to the developers.")
        end
        save(BinaryFuseFilterDBperChrom(dbi, bff_fwd, bff_rve, chrom_name), 
             file_path(storage_path, "BinaryFuseFilterDB_" * chrom_name * ".bin"))
    end

    # TODO add ambiguity handling here?!
    ambig = nothing # ambig = length(ambig) > 0 ? AmbigIdx(ambig, nothing) : nothing
    close(ref)

    save(BinaryFuseFilterDB(dbi, mpt, ambig), 
         file_path(storage_path, "BinaryFuseFilterDB.bin"))

    @info "Finished."
    return 
end




function search_chrom2(
    chrom::String, 
    detail::String, 
    guides::Vector{LongDNA{4}},
    bffddbir::String, 
    fmidbdir::String, 
    bffDB::BinaryFuseFilterDB,
    distance::Int)

    guides_uint64 = guide_to_template_format.(copy(guides))
    guides_uint64_rc = guide_to_template_format.(reverse_complement.(copy(guides)))
    guides_fmi = guide_to_template_format.(copy(guides_); alphabet = ALPHABET_UINT8)
    guides_fmi_rc = guide_to_template_format.(reverse_complement.(copy(guides_)); alphabet = ALPHABET_UINT8)


    #ref = open(genomepath, "r")
    #reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    #seq = getchromseq(gi.is_fa, reader[chrom])
    fmi = load(joinpath(storage_dir, chrom * ".bin"))
    detail_path = joinpath(detail, "detail_" * string(chrom) * ".csv")
    detail_file = open(detail_path, "w")

    fwd_offt = Vector{Vector{Path}}()
    for (i, s) in enumerate(guides)

        matp = guide_in_template_format[pathTemplate]
        map(as64, eachrow(matp))
        
        # TODO - we can chromosome specific hashes - after all each chromosome is on its own when we perform the search
        pat = ARTEMIS.templates_to_sequences(s, mpt; dist = distance)

        pat_in_genome = falses(lastindex(pat))
        for (o, ot) in enumerate(pat)
            subs = expand_ambiguous(LongDNA{4}(ot.seq) * repeat(dna"N", len - length(ot.seq)))
            for sub in subs 
                if ARTEMIS.get_count_idx(db.bins, convert(UInt64, sub), right) == 0
                    pat_in_genome[o] = true
                    continue # we skip the checks as one of the subsequences was found in the genome
                end
            end
        end

        @info "Fraction of seqeunces that passed " * string(ceil(sum(pat_in_genome)/lastindex(pat_in_genome); digits = 4))
        push!(fwd_offt, pat[pat_in_genome]) # TODO probably replace with preallocated 
    end


    # wroking on this guide and his all possible off-targets
    for (i, fwd_offt_i) in enumerate(fwd_offt) # for each guide

        for offt in fwd_offt_i # for each potential OT

            fwd_iter = locate(offt.seq, fmi)
            if motif.extends5 
                fwd_iter = fwd_iter .+ length(offt.seq) .- 1
            end
            for pos in fwd_iter # TODO this might be slow too
                line = string(guides[i]) * "," * "no_alignment" * "," * 
                    string(offt.seq) * "," * string(offt.dist) * "," *
                    chrom * "," * string(pos) * "," * "+" * "\n"
                write(detail_file, line)
            end
                
            rve_iter = locate(reverse_complement(offt.seq), fmi)
            if !motif.extends5
                rve_iter = rve_iter .+ length(offt.seq) .- 1
            end
            for pos in rve_iter
                line = string(guides[i]) * "," * "no_alignment" * "," * 
                    string(offt.seq) * "," * string(offt.dist) * "," *
                    chrom * "," * string(pos) * "," * "-" * "\n"
                write(detail_file, line)
            end
        end
    end
    close(detail_file)
    #close(ref)
    

    
    return 
end



"""
```
search_BinaryFuseFilterDB(
    bffddbir::String, 
    fmidbdir::String,
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 4)
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
function search_BinaryFuseFilterDB(
    bffddbir::String, 
    fmidbdir::String,
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3)

    if any(isambig.(guides)) # TODO?
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    guides_ = copy(guides)
    len_noPAM_noEXT = length_noPAM(motif)
    len = len_noPAM_noEXT + motif.distance

    if any(len_noPAM_noEXT .!= length.(guides_))
        throw("Guide queries are not of the correct length to use with this Motif: " * string(db.dbi.motif))
    end

    # reverse guides so that PAM is always on the left
    #if db.dbi.motif.extends5
    #    guides_ = reverse.(guides_)
    #end

    bffDB = load(joinpath(bffddbir, "binaryFuseFilterDB.bin"))
    bffDB = restrictDistance!(bffDB, distance)
    # we input guides that are in forward search configuration # which means 20bp-NGG
    ThreadsX.map(ch -> search_chrom2(ch, dirname(output_file), guides_, bffddbir, fmidbdir, bffDB, distance), bffDB.dbi.gi.chrom)
    cleanup_detail(output_file)
    return 
end