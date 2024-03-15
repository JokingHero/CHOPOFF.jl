struct SymbolicAlignments # PathTemplates, but without the fluff and length restricted 
    dbi::DBInfo # has Motif
    paths::Matrix{Int}
    paths_distances::Vector{Int}
    hash_len::Int
end

struct AmbigPrefixHashDB
    mpt::Union{Nothing, SymbolicAlignments} # Nothing when Motif allows ambig in PrefixHashDB

    prefix_hash::Vector{<:Unsigned} # sorted, hashes can repeat and point to different groupings
    prefix_start::Vector{<:Unsigned} # grouping pointers - start + width (smaller)
    prefix_width::Vector{<:Unsigned}

    suffix::Vector{<:Unsigned} # leftover of the hash, as small as possible
    chrom::Vector{<:Unsigned}
    pos::Vector{<:Unsigned}
    isplus::Vector{Bool}
    annot::Union{Vector{AbstractString}, Nothing} # use InlineStrings here to save time loading it up
    is_ambig::BitVector # because length of the OTs is known we can have each position showing ambiguity value
    ambig::Vector{UInt8} # encoding of special bases
end

struct PrefixHashDB
    mpt::SymbolicAlignments
    
    prefix::Vector{<:Unsigned} # sorted, unique 
    #prefix_start::Vector{<:Unsigned} # pointers to suffix
    #prefix_width::Vector{<:Unsigned}

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
    ambig_guides::AG,
    ambig_chrom::C,
    ambig_pos::P,
    ambig_isplus::Vector{Bool},
    guides::H,
    chrom::C,
    pos::P,
    isplus::Vector{Bool},
    chrom_name::CU,
    dbi::DBInfo,
    chrom_seq::K,   
    is_antisense::Bool) where {
        AG<:Vector{LongDNA{4}}, C<:Vector{<:Unsigned}, P<:Vector{<:Unsigned}, 
        H<:Vector{<:Unsigned}, CU<:Unsigned, K<:BioSequence}

    pam_loci = is_antisense ? dbi.motif.pam_loci_rve : dbi.motif.pam_loci_fwd

    if length(dbi.motif) != 0
        this_pos = ARTEMIS.findguides(dbi, chrom_seq, is_antisense)
        this_guides = ThreadsX.map(x -> ARTEMIS.removepam(chrom_seq[x], pam_loci), this_pos)
        if dbi.motif.distance > 0
            this_guides = ARTEMIS.add_extension(this_guides, this_pos, dbi, chrom_seq, is_antisense)
        end
        this_guides, this_pos = ARTEMIS.normalize_to_PAMseqEXT(this_guides, this_pos, dbi, is_antisense)
        this_pos = ARTEMIS.convert.(dbi.gi.pos_type, this_pos)

        idx = ThreadsX.map(ARTEMIS.isambig, this_guides)
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
```julia-repl
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
    hash_type = ARTEMIS.smallestutype(parse(UInt, repeat("1", hash_len * 2); base = 2))
    suffix_type = ARTEMIS.smallestutype(parse(UInt, repeat("1", (ARTEMIS.length_noPAM(motif) - hash_len + motif.distance) * 2); base = 2))
    ot_type = ARTEMIS.smallestutype(parse(UInt, repeat("1", (ARTEMIS.length_noPAM(motif) + motif.distance) * 2); base = 2))

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

    # step 2
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

    if length(ambig_guides) > 0 # TODO
        @info "Step 4: Constructing DB for ambigous gRNAs."
        #=
        ambig_guides
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
        =#
    end
    
    return storage_dir
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

Find all off-targets for `guides` within distance of `dist` using linearHashDB located at `storage_dir`.
Uses early stopping to stop searching when a guide passes a limit on number of off-targets. This method does not 
keep track of the off-target locations and does not filter overlapping off-targets, therefore it might hit the 
early stopping condition a little earlier than intended.

Assumes your guides do not contain PAM, and are all in the same direction as 
you would order from the lab e.g.:

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
```julia-repl
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
    ot_len = ARTEMIS.length_noPAM(db.mpt.dbi.motif) + db.mpt.dbi.motif.distance
    ot_type = ARTEMIS.smallestutype(parse(UInt, repeat("1", ot_len * 2); base = 2))
    s_len = ot_len - db.mpt.hash_len

    if distance > db.mpt.dbi.motif.distance
        error("For this database maximum distance is " * string(db.mpt.dbi.motif.distance))
    end

    if db.mpt.dbi.motif.ambig_max > 0
        # TODO
    end

    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if db.mpt.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    paths = db.mpt.paths[db.mpt.paths_distances .<= distance, :]
    mkpath(dirname(output_file))

    ThreadsX.map(guides_) do g
        guides_formated = ARTEMIS.guide_to_template_format(g; alphabet = ARTEMIS.ALPHABET_TWOBIT)
        sa = guides_formated[paths]
        sa = Base.map(x -> ARTEMIS.asUInt(eltype(db.prefix), x), eachrow(sa))
        sa = Set(sa)
        sa = in.(db.prefix, Ref(sa))
        
        es_acc = zeros(Int64, length(early_stopping))
        is_es = false
        ots = LongDNA{4}.((convert.(ot_type, db.prefix[sa]) .<< (2 * s_len)) .| convert.(ot_type, db.suffix[sa]), ot_len)
        
        detail_path = joinpath(dirname(output_file), "detail_" * string(g) * ".csv")
        detail_file = open(detail_path, "w")
        guide_stranded = db.mpt.dbi.motif.extends5 ? reverse(g) : g
        guide_stranded = string(guide_stranded)
        if length(ots) == 0
            close(detail_file)
            return
        end

        aln = ARTEMIS.align(g, ots[1], distance, iscompatible)
        if aln.dist <= distance
            es_acc[aln.dist + 1] += 1
            if db.mpt.dbi.motif.extends5
                aln_guide = reverse(aln.guide)
                aln_ref = reverse(aln.ref)
            else
                aln_guide = aln.guide
                aln_ref = aln.ref
            end
            strand = db.isplus[sa][1] ? "+" : "-"
            ot = guide_stranded * "," * aln_guide * "," * 
                aln_ref * "," * string(aln.dist) * "," *
                db.mpt.dbi.gi.chrom[db.chrom[sa][1]] * "," * 
                string(db.pos[sa][1]) * "," * strand * "\n"
            write(detail_file, ot)
        end
        for i in 2:length(ots)
            ot = ots[i]
            if ot != ots[i - 1] # we can recycle alignments because ots are sorted
                aln = ARTEMIS.align(g, ot, distance, iscompatible)
            end
            
            if aln.dist <= distance
                if db.mpt.dbi.motif.extends5
                    aln_guide = reverse(aln.guide)
                    aln_ref = reverse(aln.ref)
                else
                    aln_guide = aln.guide
                    aln_ref = aln.ref
                end
                strand = db.isplus[sa][i] ? "+" : "-"
                ot = guide_stranded * "," * aln_guide * "," * 
                    aln_ref * "," * string(aln.dist) * "," *
                    db.mpt.dbi.gi.chrom[db.chrom[sa][i]] * "," * 
                    string(db.pos[sa][i]) * "," * strand * "\n"
                write(detail_file, ot)
                es_acc[aln.dist + 1] += 1
                if es_acc[aln.dist + 1] >= early_stopping[aln.dist + 1]
                    is_es = true
                    break
                end
            end
        end
        close(detail_file)
        return
    end

    cleanup_detail(output_file)
    return
end