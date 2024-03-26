struct SymbolicAlignments # PathTemplates, but without the fluff and length restricted 
    dbi::DBInfo # has Motif
    paths::Matrix{Int}
    paths_distances::Vector{Int}
    hash_len::Int
end

struct AmbigPrefixHashDB
    mpt::Union{Nothing, CHOPOFF.SymbolicAlignments} # Nothing when Motif allows ambig in PrefixHashDB
    prefix::Vector{<:Unsigned} # sorted, hashes can repeat and point to different groupings
    prefix_idx::Vector{Int}
    suffix::Vector{<:Unsigned} # leftover of the hash - ambigous bases are set to 11
    chrom::Vector{<:Unsigned}
    pos::Vector{<:Unsigned}
    isplus::Vector{Bool}
    annot::Union{Vector{<:AbstractString}, Nothing}
    is_ambig::BitVector # because length of the OTs is known we can have each position showing ambiguity value
    ambig::Vector{UInt8} # encoding of special bases
end

function build_ambigPrefixHashDB(
    ambig_guides::Vector{LongDNA{4}},
    ambig_chrom::Vector{<:Unsigned},
    ambig_pos::Vector{<:Unsigned},
    ambig_isplus::Vector{Bool},
    l::Int, # length_noPAM + distance
    hash_len::Int,
    ot_type::Type, # type of the OTs
    hash_type::Type,
    suffix_type::Type,
    annot::Union{Vector{<:AbstractString}, Nothing},
    mpt::Union{Nothing, CHOPOFF.SymbolicAlignments})

    ambig_prefixes = Vector{ot_type}()
    ambig_indexes = Vector{Int}()
    ambig_suffixes = Vector{ot_type}()
    ambig_bv = BitVector()
    ambig_bases = Vector{UInt8}()
    for (i, ag) in enumerate(ambig_guides)
        prefixes, idx = CHOPOFF.expand_ambiguous(ag) 
        prefixes = convert.(ot_type, prefixes)
        append!(ambig_suffixes, prefixes[1]) # any of the prefixes will do
        prefixes = convert.(hash_type, prefixes .>> (2 * (l - hash_len)))
        prefixes = unique(prefixes)
        append!(ambig_prefixes, prefixes)
        append!(ambig_indexes, repeat([i], length(prefixes)))
        bv = BitVector(falses(l))
        bv[idx] .= 1
        append!(ambig_bv, bv)
        append!(ambig_bases, convert.(UInt8, ag[idx]))
    end

    ambig_guides = nothing
    mask = (one(ot_type) << (2 * (l - hash_len))) - one(ot_type)
    ambig_suffixes = convert.(suffix_type, mask .& ambig_suffixes)
    order = sortperm(ambig_prefixes)
    ambig_prefixes = ambig_prefixes[order]
    ambig_indexes = ambig_indexes[order]
    ambig_isplus = BitVector(ambig_isplus)

    return AmbigPrefixHashDB(
        mpt,
        ambig_prefixes,
        ambig_indexes,
        ambig_suffixes,
        ambig_chrom,
        ambig_pos,
        ambig_isplus,
        annot,
        ambig_bv,
        ambig_bases)
end

struct PrefixHashDB
    mpt::SymbolicAlignments
    prefix::Vector{<:Unsigned} # sorted, unique 
    suffix::Vector{<:Unsigned}
    chrom::Vector{<:Unsigned}
    pos::Vector{<:Unsigned}
    isplus::Vector{Bool}
end

struct OffTargets
    guide::LongDNA{4}
    is_es::Bool
    alignments::Vector{Aln}
    chrom::Vector{<:Unsigned}
    pos::Vector{<:Unsigned}
    isplus::Vector{Bool}
end

function gather_hashes!(
    ambig_guides::Vector{LongDNA{4}},
    ambig_chrom::Vector{<:Unsigned},
    ambig_pos::Vector{<:Unsigned},
    ambig_isplus::Vector{Bool},
    guides::Vector{<:Unsigned},
    chrom::Vector{<:Unsigned},
    pos::Vector{<:Unsigned},
    isplus::Vector{Bool},
    chrom_name::CU,
    dbi::DBInfo,
    chrom_seq::K,   
    is_antisense::Bool) where {CU<:Unsigned, K<:BioSequence}

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd

    if length(dbi.motif) != 0
        this_pos = CHOPOFF.findguides(dbi, chrom_seq, is_antisense)
        this_guides = ThreadsX.map(x -> CHOPOFF.removepam(chrom_seq[x], pam_loci), this_pos)
        if dbi.motif.distance > 0
            this_guides = CHOPOFF.add_extension(this_guides, this_pos, dbi, chrom_seq, is_antisense)
        end
        this_guides, this_pos = CHOPOFF.normalize_to_PAMseqEXT(this_guides, this_pos, dbi, is_antisense)
        this_pos = CHOPOFF.convert.(dbi.gi.pos_type, this_pos)

        idx = ThreadsX.map(CHOPOFF.isambig, this_guides)
        idx_sum = sum(idx)
        if sum(idx) != 0
            append!(ambig_guides, this_guides[idx])
            this_guides = this_guides[.!idx]
            append!(ambig_chrom, repeat([chrom_name], idx_sum))
            append!(ambig_pos, this_pos[idx])
            this_pos = this_pos[.!idx]
            append!(ambig_isplus, repeat([!is_antisense], idx_sum))
        end

        this_guides = convert.(eltype(guides), this_guides)
        append!(guides, this_guides)
        append!(chrom, repeat([chrom_name], length(this_pos)))
        append!(pos, this_pos)
        append!(isplus, repeat([!is_antisense], length(this_pos)))
    end
    return nothing
end


"""
```
name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    hash_len::Int = min(length_noPAM(motif) - motif.distance, 16))
```

Prepare prefixHashDB index for future searches using `search_prefixHashDB`.

Will return a path to the database location, same as `storage_dir`.
If interested with searches within distance 4, preferably use `prefix_len` of 8 or 9.
You can also play with `hash_len` parameter, but keeping it at 16 should be close to optimal.

# Arguments
`name` - Your preferred name for this index for easier identification.

`genomepath` - Path to the genome file, it can either be fasta or 2bit file. In case of fasta
               also prepare fasta index file with ".fai" extension.

`motif`   - Motif defines what kind of gRNA to search for.

`storage_dir`  - Folder path to the where index will be saved with name `linearDB.bin` and many prefix files.

`hash_len` - Length of the hash in bp. At maximum 16.

# Examples
```julia
$(make_example_doc("prefixHashDB"))
```
"""
function build_prefixHashDB(
    name::String,
    genomepath::String,
    motif::Motif,
    storage_dir::String,
    hash_len::Int = min(length_noPAM(motif) - motif.distance, 16))

    if hash_len > 16
        throw("hash_len $hash_len is more than 16")
    end

    dbi = DBInfo(genomepath, name, motif)
    hash_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    suffix_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", (CHOPOFF.length_noPAM(motif) - hash_len + motif.distance) * 2); base = 2))
    ot_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", (CHOPOFF.length_noPAM(motif) + motif.distance) * 2); base = 2))

    # step 1
    @info "Step 1: Searching chromosomes."
    ref = open(dbi.gi.filepath, "r")
    reader = dbi.gi.is_fa ? FASTA.Reader(ref, index = dbi.gi.filepath * ".fai") : TwoBit.Reader(ref)

    ambig_guides = Vector{LongDNA{4}}()
    ambig_chrom = Vector{dbi.gi.chrom_type}()
    ambig_pos = Vector{dbi.gi.pos_type}()
    ambig_isplus = Vector{Bool}()
    guides = Vector{ot_type}()
    chrom = Vector{dbi.gi.chrom_type}()
    pos = Vector{dbi.gi.pos_type}()
    isplus = Vector{Bool}()

    for chrom_name in dbi.gi.chrom
        record = reader[chrom_name]
        chrom_name_ = convert(dbi.gi.chrom_type, findfirst(isequal(chrom_name), dbi.gi.chrom))
        @info "Working on $chrom_name"
        chrom_seq = dbi.gi.is_fa ? FASTA.sequence(LongDNA{4}, record) : TwoBit.sequence(LongDNA{4}, record)
        gather_hashes!(
            ambig_guides, ambig_chrom, ambig_pos, ambig_isplus,
            guides, chrom, pos, isplus, 
            chrom_name_, dbi, chrom_seq, false)
        gather_hashes!(
            ambig_guides, ambig_chrom, ambig_pos, ambig_isplus,
            guides, chrom, pos, isplus, 
            chrom_name_, dbi, chrom_seq, true)
    end

    @info "Step 2: Constructing DB."
    # split gRNA into prefix - suffix
    l = length_noPAM(motif) + motif.distance
    prefixes = convert.(hash_type, guides .>> (2 * (l - hash_len)))
    mask = (one(ot_type) << (2 * (l - hash_len))) - one(ot_type)
    suffixes = convert.(suffix_type, mask .& guides)
    guides = nothing

    order = sortperm(prefixes)
    prefixes = prefixes[order]
    suffixes = suffixes[order]
    chrom = chrom[order]
    pos = pos[order]
    isplus = BitVector(isplus[order])

    @info "Step 3: Constructing Paths for hashes"
    mpt = build_PathTemplates(motif; restrict_to_len = hash_len, withPAM = false)
    paths = mpt.paths[:, 1:hash_len]
    not_dups = map(!, BitVector(nonunique(DataFrame(paths, :auto)))) # how can there be no duplicated function?!
    mkpath(storage_dir)
    save(PrefixHashDB(SymbolicAlignments(dbi, paths[not_dups, :], mpt.distances[not_dups], hash_len),
        prefixes, suffixes, chrom, pos, isplus), joinpath(storage_dir, "prefixHashDB.bin"))
    @info "Finished constructing prefixHashDB in " * storage_dir

    if length(ambig_guides) > 0
        @info "Step 4: Constructing DB for ambigous gRNAs."
        paths = nothing
        mpt = nothing
        prefixes = nothing
        suffixes = nothing
        chrom = nothing
        pos = nothing
        isplus = nothing
        
        save(build_ambigPrefixHashDB(
                ambig_guides, ambig_chrom, ambig_pos, ambig_isplus,
                l, hash_len, ot_type, hash_type, suffix_type, nothing, nothing), 
            joinpath(storage_dir, "ambigPrefixHashDB.bin"))
        @info "Finished constructing DB for ambigous OTs in " * storage_dir
    end
    
    return storage_dir
end


# sa - sometimes small, sometimes large
# prefixes - around 300M                       
function potential_ots_idx(sa::Vector{<:Unsigned}, prefixes::Vector{<:Unsigned})
    idx = searchsorted.(Ref(prefixes), sa)
    return filter(x -> x.start <= x.stop, idx)
end

"""
```
search_prefixHashDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    early_stopping::Vector{Int} = Int.(floor.(exp.(0:distance))))
```

Find all off-targets for `guides` within distance of `dist` using prefixHashDB located at `storage_dir`.
Uses early stopping to stop searching when a guide passes a limit on number of off-targets. This method does not 
keep track of the off-target locations and does not filter overlapping off-targets, therefore it might hit the 
early stopping condition a little earlier than intended.

Assumes your guides do not contain PAM, and are all in the same direction as 
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
$(make_example_doc("prefixHashDB"))
```
"""
function search_prefixHashDB(
    storage_dir::String, 
    guides::Vector{LongDNA{4}}, 
    output_file::String;
    distance::Int = 3,
    early_stopping::Vector{Int} = Int.(floor.(exp.(0:distance))))

    if length(early_stopping) != (distance + 1)
        error("Specify one early stopping condition for a each distance, starting from distance 0.")
    end

    db = load(joinpath(storage_dir, "prefixHashDB.bin"))
    ot_len = CHOPOFF.length_noPAM(db.mpt.dbi.motif) + db.mpt.dbi.motif.distance
    ot_type = CHOPOFF.smallestutype(parse(UInt, repeat("1", ot_len * 2); base = 2))
    s_len = ot_len - db.mpt.hash_len

    if distance > db.mpt.dbi.motif.distance
        error("For this database maximum distance is " * string(db.mpt.dbi.motif.distance))
    end

    adb_file = joinpath(storage_dir, "ambigPrefixHashDB.bin")
    if isfile(adb_file)
        adb = load(adb_file)
    else 
        adb = nothing
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if db.mpt.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    paths = db.mpt.paths[db.mpt.paths_distances .<= distance, :]
    mkpath(dirname(output_file))

    Base.map(guides_) do g # maybe a function would be faster than lambda here?
        guides_formated = CHOPOFF.guide_to_template_format(g; alphabet = CHOPOFF.ALPHABET_TWOBIT)
        sa = guides_formated[paths]
        sa = Base.map(x -> CHOPOFF.asUInt(eltype(db.prefix), x), eachrow(sa))
        sa = unique(sa)
        if !isnothing(adb)
            asa = CHOPOFF.potential_ots_idx(sa, adb.prefix)
        else
            asa = []
        end
        sa = CHOPOFF.potential_ots_idx(sa, db.prefix)
        
        es_acc = zeros(Int64, length(early_stopping))   
        detail_path = joinpath(dirname(output_file), "detail_" * string(g) * ".csv")
        detail_file = open(detail_path, "w")
        guide_stranded = db.mpt.dbi.motif.extends5 ? reverse(g) : g
        guide_stranded = string(guide_stranded)
        if length(sa) == 0 && length(asa) == 0
            close(detail_file)
            return
        end

        @inbounds for i in sa # each sa is range of indices of prefixes where all ots are the same
            ot = LongDNA{4}((convert(ot_type, db.prefix[i.start]) << (2 * s_len)) | 
                convert(ot_type, db.suffix[i.start]), ot_len)
            aln = CHOPOFF.align(g, ot, distance, iscompatible)
            if aln.dist <= distance
                if db.mpt.dbi.motif.extends5
                    aln_guide = reverse(aln.guide)
                    aln_ref = reverse(aln.ref)
                else
                    aln_guide = aln.guide
                    aln_ref = aln.ref
                end
                @inbounds for idx in i
                    strand = db.isplus[idx] ? "+" : "-"
                    ot = guide_stranded * "," * aln_guide * "," * 
                        aln_ref * "," * string(aln.dist) * "," *
                        db.mpt.dbi.gi.chrom[db.chrom[idx]] * "," * 
                        string(db.pos[idx]) * "," * strand * "\n"
                    write(detail_file, ot)
                    es_acc[aln.dist + 1] += 1
                    if es_acc[aln.dist + 1] >= early_stopping[aln.dist + 1]
                        close(detail_file)
                        return
                    end
                end
            end
        end

        if length(asa) > 0
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
                    if db.mpt.dbi.motif.extends5
                        aln_guide = reverse(aln.guide)
                        aln_ref = reverse(aln.ref)
                    else
                        aln_guide = aln.guide
                        aln_ref = aln.ref
                    end
                    strand = adb.isplus[idx] ? "+" : "-"
                    ot = guide_stranded * "," * aln_guide * "," * 
                        aln_ref * "," * string(aln.dist) * "," *
                        db.mpt.dbi.gi.chrom[adb.chrom[idx]] * "," * 
                        string(adb.pos[idx]) * "," * strand * "\n"
                    write(detail_file, ot)
                    es_acc[aln.dist + 1] += 1
                    if es_acc[aln.dist + 1] >= early_stopping[aln.dist + 1]
                        close(detail_file)
                        return
                    end
                end
            end
        end

        close(detail_file)
        return
    end

    cleanup_detail(output_file)
    return
end