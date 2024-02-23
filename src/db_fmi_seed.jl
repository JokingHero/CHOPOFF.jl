# Position should be the position of the first/last base of PAM
function make_line(guide::LongDNA{4}, chrom::String, oriented_pos::Int, aln::Aln, rev::Bool, strand::String)
    aln_guide = aln.guide
    aln_ref = aln.ref
    if rev
        aln_guide = reverse(aln_guide)
        aln_ref = reverse(aln_ref)
    end
    return string(guide) * "," * aln_guide * "," * 
        aln_ref * "," * string(aln.dist) * "," *
        chrom * "," * string(oriented_pos) * "," * strand * "\n"
end


#          NGG  CCN  TTTN NAAA
# align 
# rev       T   F    F     T
# compl     F   T    F     T
#          GGN  GGN  TTTN  TTTTN

# print         
# rev       T   T    F     F
# compl     F   F    F     F
function align_and_write!(guide::LongDNA{4}, ref::LongDNA{4}, chrom::String, start::Int, stop::Int, 
    is_antisense::Bool, motif::Motif, filter::Set{Int64}, detail_file::IOStream)
    # guide here is in program input state
    g_PAM = appendPAM_forward(guide, motif)
    # ref here is sense/antisense specific
    pam_oriented_pos = start
    if motif.extends5 && is_antisense
        # CCN ... EXT - ref
        complement!(ref)
        # becomes GGN ... EXT - ref
        reverse!(g_PAM)
    elseif motif.extends5 && !is_antisense 
        # EXT ... NGG - ref
        reverse!(ref)
        # becomes GGN ... EXT
        reverse!(g_PAM)
        pam_oriented_pos = stop
    else# !motif.extends5 && is_antisense 
        # EXT ... NAA
        reverse_complement!(ref)
        # becomes TTA ... EXT
        pam_oriented_pos = stop
        #!motif.extends5 && !is_antisense
        # TTN ... EXT
    end

    # we need PAM-guide for the alignment
    aln = align(g_PAM, ref, motif.distance)
    if aln.dist <= motif.distance
        push!(filter, pam_oriented_pos)
        
        write(detail_file, make_line(
            guide, 
            chrom,
            pam_oriented_pos,
            aln,
            motif.extends5,
            is_antisense ? "-" : "+"))
            
    end
end

# is_start indicates whether pos is anchored on the left or on the right
# simplest case possible here
function pam00!(
    guide::LongDNA{4}, seq::LongDNA{4}, chrom::String, pos::Int, is_start::Bool, is_antisense::Bool, 
    motif::Motif, filter::Set{Int64}, detail_file::IOStream)
    if !(pos in filter)
        if is_start
            start = pos
            stop = start + length(motif) + motif.distance - 1
        else
            start = pos - length(motif) - motif.distance - 1
            stop = pos
        end
        align_and_write!(guide, seq[start:stop], chrom, start, stop, is_antisense, motif, filter, detail_file)
    end
end


function pam100_anchor_right!(
    guide::LongDNA{4}, seq::LongDNA{4}, chrom::String,
    pos::Int, pam_len::Int, err_drift::Int,
    pam_pos::Vector{Int}, 
    is_antisense::Bool, motif::Motif, 
    filter::Set{Int64}, detail_file::IOStream)
    # find first value >= to p
        # worst case 00111PAM
        #            ^
        #            10011PAM
        #             ^
        # best case  11001PAM
        #              ^
        # we search with -err_drift to get first possible legal PAM on the maximal left err_drift
        stop_nopam = pos - pam_len - 1 - err_drift
        idx = searchsortedfirst(pam_pos, stop_nopam)
    
        # 00111-30+3PAM
        #     ^     ^
        while idx <= length(pam_pos) && (pam_pos[idx] - stop_nopam <= (err_drift * 2 + 1))
            # we have to check this PAM idx
            pam00!(guide, seq, chrom, pam_pos[idx] + pam_len - 1, false, is_antisense, 
                motif, filter, detail_file)
            # now that we know this PAM passed the proximity test
            # this is not necessarily optimal PAM for this search
            # we need to check upstream for potentially more convenient PAMs...
            idx += 1
        end
end


function pam100_anchor_left!(
    guide::LongDNA{4}, seq::LongDNA{4}, chrom::String,
    pos::Int, pam_len::Int, err_drift::Int, pam_pos::Vector{Int}, 
    is_antisense::Bool, motif::Motif, 
    filter::Set{Int64}, detail_file::IOStream)
    # worst case PAM11100
    #                  ^
    #            PAM11001
    #                 ^
    #            PAM10011
    #                ^
    #            12345678
    start_nopam = pos + pam_len + err_drift
    # find last value <= to p
    idx = searchsortedlast(pam_pos, start_nopam)
    # PAM+30-311100
    #   ^        ^
    while idx > 0 && idx <= length(pam_pos) && 
        (start_nopam - pam_pos[idx] <= (err_drift * 2 + 1))
        # we have to check this PAM idx

        pam00!(guide, seq, chrom, pam_pos[idx] - pam_len + 1, true, 
            is_antisense, motif, filter, detail_file)
        idx -= 1
    end
end


"""
```
search_fmiDB_seed(
    guides::Vector{LongDNA{4}},
    fmidbdir::String, 
    genomepath::String, 
    pamDB::PAMinFMI,
    output_file::String; 
    distance::Int = 2)
```

Search FM-index for off-targets using 01*0 seed method. 
Read more here: [publication](https://www.sciencedirect.com/science/article/pii/S1570866716300028) 
and [pdf](https://mikael-salson.univ-lille.fr//articles/VST_Iwoca14.pdf).

**Experimental! Proof-of-concept!**


# Arguments
`guides` - a vector of gRNAs without PAM.

`fmidbdir`   - Path to the folder where FM-index was build using `build_fmi`.

`genomepath` - Path to the genome used to build the FM-index.

`pamDB` - object build with `build_pamDB` that contains locations of the PAM inside the genome.

`output_file`  - Where output will be saved.

`distance`  - Search distance, maximum of 2 is practical.


# Examples
```julia-repl
# prepare libs
using ARTEMIS, BioSequences

# make a temporary directory
tdir = tempname()
fmi_dir = joinpath(tdir, "fmi")
mkpath(fmi_dir)

# use ARTEMIS example genome
artemis_path = splitpath(dirname(pathof(ARTEMIS)))[1:end-1]
genome = joinpath(vcat(artemis_path, 
    "test", "sample_data", "genome", "semirandom.fa"))
# build FM-index
build_fmiDB(genome, fmi_dir)

# build a pamDB
motif = Motif("Cas9"; distance = 1)
pamDB = build_pamDB(fmi_dir, motif)

# prepare PathTemplates
mpt = build_PathTemplates(motif)

# prepare output folder
res_dir = joinpath(tdir, "results")
mkpath(res_dir)

# load up example gRNAs
guides_s = Set(readlines(joinpath(vcat(artemis_path, 
    "test", "sample_data", "crispritz_results", "guides.txt"))))
guides = LongDNA{4}.(guides_s)
    
# finally, make results!
res_path = joinpath(res_dir, "results.csv")
search_fmiDB_seed(guides, fmi_dir, genome, pamDB, res_path; distance = 1)

# load results
using DataFrames, CSV
res = DataFrame(CSV.File(res_path))

# filter results by close proximity
res = filter_overlapping(res, 23)

# summarize results into a table of counts by distance
summary = summarize_offtargets(res; distance = 1)
```
"""
function search_fmiDB_seed(guides::Vector{LongDNA{4}},
    fmidbdir::String, genomepath::String, pamDB::PAMinFMI,
    output_file::String; distance::Int = 2)

    motif = setdist(pamDB.motif, distance)
    detail_path = dirname(output_file)
    chunks = distance + 2
    combs = collect(combinations(1:chunks, 2))
    # reorder combs so that options with PAM00 are first, then PAM0 then other
    if motif.extends5
        combs_rev = combs
        combs = sort(combs, rev = true, by = x -> (x[2], x[1]))
    else
        throw("Sorry, but not ready for testing of extends5 false yet.")
        combs_rev = sort(combs, rev = true, by = x -> (x[2], x[1]))
    end

    gi2 = load(joinpath(fmidbdir, "genomeInfo.bin"))
    gi = GenomeInfo(genomepath)
    if !isequal(gi, gi2)
        msg = "Supplied genome is different than the genome used for building of pamDB. "
        if isequal(Set(gi.chrom), Set(gi.chrom))
            msg *= "Supplied genome has the same chromosome names."
        else
            if all(occursin.(gi.chrom, gi2.chrom))
                msg *= "Supplied genome has different chromosome names, but all exist in pamDB."
            else
                throw("Supplied genome has different chromosome names, but not all exist in pamDB.")
            end
        end
        @warn msg
    end

    guides = copy(guides)
    chunk_len = Int(floor(length(guides[1]) / chunks)) # this can not be used in actual calculations because chunks might not be of equal size...

    # important here is that last chunk is dependant on the motif Cas9 vs Cpf1 style
    # add PAM in front or in the back
    # expand ambiguous
    if motif.extends5 # Cas9 style with pattern-PAM
        l = length(guides[1]) - chunks * chunk_len
        chunk_idx = [(chunk_len*j+l+1):(chunk_len*(j+1)+l) for j = 0:chunks-1]
        chunk_idx[1]= 1:chunk_idx[1].stop # leftovers are added to the first chunk
        guide_chunks = map(x -> [x[i] for i in chunk_idx], guides)
        map(x -> x[end] = x[end] * motif.fwd[motif.pam_loci_fwd], guide_chunks)
    else # Cpf1 style with PAM-pattern
        chunk_idx = [(chunk_len*j+1):(chunk_len*(j+1)) for j = 0:chunks-1]
        chunk_idx[end]= chunk_idx[end].start:length(guides[1]) # leftovers are added to the last chunk
        guide_chunks = map(x -> [x[i] for i in chunk_idx], guides)
        # first chunk has additionally PAM in front appended
        map(x -> x[1] = motif.fwd[motif.pam_loci_fwd] * x[1], guide_chunks)
    end
    # allows for ambiguity in the guide, solves ambiguity in the PAM
    guide_chunks = map(guide_chunks) do gch
        map(gch) do x
            expand_ambiguous(x)
        end
    end
    pam_len = length(motif) - length_noPAM(motif)

    # we input guides that are in forward search configuration
    ref = open(genomepath, "r")
    reader = gi.is_fa ? FASTA.Reader(ref, index=genomepath * ".fai") : TwoBit.Reader(ref)
    if Threads.nthreads() != 1
        reader = collect(reader) # this is potentially explosive in terms of memory
        # TODO parallelization without pre-loading each reader, maybe with ThreadsX.mapi
    end

    ThreadsX.foreach(enumerate(gi.chrom)) do (ic, chrom)
        seq = getchromseq(gi.is_fa, Threads.nthreads() == 1 ? reader[chrom] : reader[ic])
        fmi = load(joinpath(fmidbdir, chrom * ".bin"))
        pam_fwd = pamDB.pam_loc_fwd[ic]
        pam_rve = pamDB.pam_loc_rve[ic]

        detail_chrom = joinpath(detail_path, "detail_" * string(chrom) * ".csv")
        detail_file = open(detail_chrom, "w")

        # working on this guide and his all possible off-targets
        for (i, g_chunks) in enumerate(guide_chunks) # for each guide
            g_chunks_len = map(x -> length(x[1]), g_chunks)
            g_chunks_len_rev = reverse(g_chunks_len)
            g_chunks_rev = reverse(map(x -> reverse_complement.(x), g_chunks))
            
            # we want to have PAM searches first because this will guarantee finding 0-distance first
            # FORWARD
            # keep track of all found so far to not search same PAM multiple times
            fwd_filter = Set{Int64}()
            for ch in combs
                # if pam is included in the chunk, skip PAM matching
                pam_in_chunk = (motif.extends5 && ch[2] == length(g_chunks)) || (!motif.extends5 && ch[1] == 1)
                # if two chunks are adjacent - merge and skip distance matching
                single_ch = (ch[2] - ch[1]) == 1
                if single_ch
                    #    00111PAM
                    #    10011PAM
                    #    11001PAM
                    #    11100PAM
                    # PAM00111
                    # PAM10011
                    # PAM11001
                    # PAM11100
                    both = map(x -> x[1] * x[2], 
                        vec(collect(Iterators.product(g_chunks[ch[1]], g_chunks[ch[2]]))))

                    # locateall returns start position of the search sequence
                    # start = end - length + 1
                    # end = start + length - 1
                    # we return as pos PAM first/last base as location in all algorithms
                    fwd_pos = mapreduce(x -> locateall(x, fmi), vcat, both)
                    both_len = length(both[1])
                    if pam_in_chunk
                        for p in fwd_pos
                            if motif.extends5
                                # on forward strand Cas9 returns last idx with PAM
                                pam00!(guides[i], seq, chrom, p + both_len - 1, false, false, motif, fwd_filter, detail_file)
                            else
                                pam00!(guides[i], seq, chrom, p, true, false, motif, fwd_filter, detail_file)
                            end
                        end
                    else # pam not in chunk..., but still two chunks are connected
                        # PAMs positions in pamDB are these:
                        #     NGG  CCN TTTN  NAAA
                        #     ^      ^    ^  ^
                        # this many -1, 0, +1 err_drifts
                        err_drift = motif.extends5 ? (length(g_chunks) - ch[2]) : (ch[1]  - 1)
                        for p in fwd_pos
                            if motif.extends5 # NGG 
                                # find first value >= to p
                                # worst case 00111PAM
                                #            ^
                                #            10011PAM
                                #             ^
                                # best case  11001PAM
                                #              ^
                                pam100_anchor_right!(guides[i], seq, chrom,
                                    p + both_len + sum(g_chunks_len[(ch[2]+1):end]), 
                                    pam_len, err_drift, pam_fwd, false, motif, fwd_filter, detail_file)
                            else # TTTN
                                pam100_anchor_left!(guides[i], seq, chrom,
                                    p - sum(g_chunks_len[1:(ch[1]-1)]) + 1, 
                                    pam_len, err_drift, pam_fwd, false, motif, fwd_filter, detail_file)
                            end
                        end
                    end
                else # not a single chunk
                    if pam_in_chunk # only the distance between chunks has to be fulfilled
                        # 01110PAM
                        # ^   ^
                        #     p
                        # 10110PAM
                        # 11010PAM
                        if motif.extends5 # we branch here because we want to iterate over chunk with PAM
                            fwd_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[1]])
                            sort!(fwd_pos_l)
                            fwd_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[2]])
                            err_drift = ch[2] - ch[1] - 1
                            for p in fwd_pos_r
                                # find closest to the left
                                stop = p + g_chunks_len[ch[2]] - 1
                                start_r = p - sum(g_chunks_len[(ch[1] + 1):(ch[2] - 1)]) + 1 + err_drift
                                idx = searchsortedlast(fwd_pos_l, start_r)
                                if idx > 0 && length(fwd_pos_l) <= idx &&
                                    ((start_r - fwd_pos_l[idx]) <= (err_drift * 2 + 1)) && !(stop in fwd_filter)
                                    start = stop - length(motif) - motif.distance + 1
                                    align_and_write!(guides[i], seq[start:stop], chrom, start, stop, false, motif, fwd_filter, detail_file)
                                end
                            end
                        else # TTTN + pam_in_chunk
                            # PAM01110
                            # ^      ^
                            # PAM01101
                            # PAM01011
                            fwd_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[1]])
                            fwd_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[2]])
                            sort!(fwd_pos_r)
                            err_drift = ch[2] - ch[1] - 1
                            for p in fwd_pos_l
                                # find closest to the left
                                start = p
                                stop_r = p + sum(g_chunks_len[1:(ch[2] - 1)]) - 1 - err_drift
                                idx = searchsortedfirst(fwd_pos_r, stop_r)
                                if idx <= length(fwd_pos_r) && 
                                    ((fwd_pos_r[idx] - stop_r) <= (err_drift * 2 + 1)) && 
                                    !(start in fwd_filter)
                                    
                                    stop = start + length(motif) + motif.distance - 1
                                    align_and_write!(guides[i], seq[start:stop], chrom, start, stop, false, motif, fwd_filter, detail_file)
                                end
                            end
                        end
                    else # distance between chunks and distance to the PAM
                        # PAMs positions in pamDB are these:
                        #     NGG  CCN TTTN  NAAA
                        #     ^      ^    ^  ^
                        
                        # 01101PAM
                        # ^  ^ ^ 
                        # 10101PAM
                        # 01011PAM
                        if motif.extends5 # we branch here because we want to iterate chunk closer to PAM first
                            fwd_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[1]])
                            sort!(fwd_pos_l)
                            fwd_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[2]])
                            err_drift = ch[2] - ch[1] - 1
                            for p in fwd_pos_r
                                # find closest chunk to the left
                                start_r = p - sum(g_chunks_len[(ch[1] + 1):(ch[2] - 1)]) + 1 + err_drift
                                idx = searchsortedlast(fwd_pos_l, start_r)
                                if idx > 0 && length(fwd_pos_l) <= idx &&
                                    ((start_r - fwd_pos_l[idx]) <= (err_drift * 2 + 1))
                                    # now we find closes upstream PAM
                                    err_pam = length(g_chunks_len) - ch[2]
                                    stop_r = p + sum(g_chunks_len[ch[2]:end]) - pam_len - 1 - err_pam
                                    # 01011-20+2PAM
                                    #     ^     ^
                                    idx = searchsortedfirst(pam_fwd, stop_r)
                                    while idx <= length(pam_fwd) && (pam_fwd[idx] - stop_r <= (err_pam * 2 + 1))
                                        pam00!(guides[i], seq, chrom, pam_fwd[idx] + pam_len - 1, false, false, motif, fwd_filter, detail_file)
                                        idx += 1
                                    end
                                end
                            end
                        else
                            # PAM10110
                            #   ^ ^  ^ 
                            # PAM10101
                            # PAM11010
                            fwd_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[1]])
                            fwd_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks[ch[2]])
                            sort!(fwd_pos_r)
                            err_drift = ch[2] - ch[1] - 1
                            for p in fwd_pos_l
                                # find closest chunk to the right
                                stop_l = p + sum(g_chunks_len[ch[1]:(ch[2]-1)]) - 1 - err_drift
                                idx = searchsortedfirst(fwd_pos_r, stop_l)
                                if (idx <= length(fwd_pos_r) && 
                                    ((fwd_pos_r[idx] - stop_l) <= (err_drift * 2 + 1)))

                                    # now find the left side PAM
                                    err_pam = ch[1] - 1
                                    start_nopam = p - sum(g_chunks_len[1:(ch[1]-1)]) + 1 + pam_len + err_pam
                                    # find last value <= to p
                                    idx = searchsortedlast(pam_fwd, start_nopam)
                                    # PAM+20-211010
                                    #   ^       ^
                                    while (idx > 0 && length(pam_fwd) <= idx 
                                        && (start_nopam - pam_fwd[idx] <= (err_pam * 2 + 1)))
                                        
                                        pam00!(guides[i], seq, chrom, pam_fwd[idx] - pam_len + 1, true, false, motif, fwd_filter, detail_file)
                                        idx -= 1
                                    end
                                end
                            end
                        end
                    end
                end
            end

            # REVERSE
            rve_filter = Set{Int64}()
            for ch in combs_rev
                # if pam is included in the chunk, skip PAM matching
                pam_in_chunk = (!motif.extends5 && ch[2] == length(g_chunks_rev)) || (motif.extends5 && ch[1] == 1)
                # if two chunks are adjacent - merge and skip distance matching
                single_ch = (ch[2] - ch[1]) == 1
                if single_ch
                    #    00111PAM
                    #    10011PAM
                    #    11001PAM
                    #    11100PAM
                    # PAM00111
                    # PAM10011
                    # PAM11001
                    # PAM11100
                    both = map(x -> x[1] * x[2], 
                        vec(collect(Iterators.product(g_chunks_rev[ch[1]], g_chunks_rev[ch[2]]))))

                    # locateall returns start position of the search sequence
                    # start = end - length + 1
                    # end = start + length - 1
                    # we return as pos PAM first/last base as location in all algorithms
                    rve_pos = mapreduce(x -> locateall(x, fmi), vcat, both)
                    both_len = length(both[1])
                    if pam_in_chunk
                        for p in rve_pos
                            if motif.extends5
                                pam00!(guides[i], seq, chrom, p, true, true, motif, rve_filter, detail_file)
                            else
                                pam00!(guides[i], seq, chrom, p + both_len - 1, false, true, motif, rve_filter, detail_file)
                            end
                        end
                    else # pam not in chunk..., but still two chunks are connected
                        # PAMs positions in pamDB are these:
                        #     NGG  CCN TTTN  NAAA
                        #     ^      ^    ^  ^
                        # this many -1, 0, +1 err_drifts
                        err_drift = motif.extends5 ? (ch[1]  - 1) : (length(g_chunks_rev) - ch[2])
                        for p in rve_pos
                            if motif.extends5
                                pam100_anchor_left!(guides[i], seq, chrom,
                                    p - sum(g_chunks_len_rev[1:(ch[1]-1)]) + 1, 
                                    pam_len, err_drift, pam_rve, true, motif, rve_filter, detail_file)
                            else
                                # worst case 00111PAM
                                #            ^
                                #            10011PAM
                                #             ^
                                # best case  11001PAM
                                #              ^
                                pam100_anchor_right!(guides[i], seq, chrom,
                                    p + both_len + sum(g_chunks_len_rev[(ch[2]+1):end]), 
                                    pam_len, err_drift, pam_rve, true, motif, rve_filter, detail_file)
                            end
                        end
                    end
                else # not a single chunk
                    if pam_in_chunk # only the distance between chunks has to be fulfilled
                        if motif.extends5 # we branch here because we want to iterate over chunk with PAM
                            # CCN + pam_in_chunk
                            # PAM01110
                            # ^      ^
                            # PAM01101
                            # PAM01011
                            rve_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[1]])
                            rve_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[2]])
                            sort!(rve_pos_r)
                            err_drift = ch[2] - ch[1] - 1
                            for p in rve_pos_l
                                # find closest to the right
                                start = p
                                stop_r = p + sum(g_chunks_len_rev[1:(ch[2] - 1)]) - 1 - err_drift
                                idx = searchsortedfirst(rve_pos_r, stop_r)
                                if idx <= length(rve_pos_r) && 
                                    ((stop_r - rve_pos_r[idx]) <= (err_drift * 2 + 1)) && 
                                    !(start in rve_filter)
                                    
                                    stop = start + length(motif) + motif.distance - 1
                                    align_and_write!(guides[i], seq[start:stop], chrom, start, stop, true, motif, rve_filter, detail_file)
                                end
                            end
                        else 
                            # 01110PAM
                            # ^   ^
                            # 10110PAM
                            # 11010PAM
                            rve_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[1]])
                            sort!(rve_pos_l)
                            rve_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[2]])
                            err_drift = ch[2] - ch[1] - 1
                            for p in rve_pos_r
                                # find closest to the left
                                stop = p + g_chunks_len_rev[ch[2]] - 1
                                start_r = p - sum(g_chunks_len_rev[(ch[1] + 1):(ch[2] - 1)]) + 1 + err_drift
                                idx = searchsortedlast(rve_pos_l, start_r)
                                if idx > 0 && length(rve_pos_l) <= idx && 
                                    ((start_r - rve_pos_l[idx]) <= (err_drift * 2 + 1)) && !(stop in fwd_filter)
                                    
                                    start = stop - length(motif) - motif.distance + 1
                                    align_and_write!(guides[i], seq[start:stop], chrom, start, stop, true, motif, rve_filter, detail_file)
                                end
                            end
                        end
                    else # distance between chunks and distance to the PAM
                        # PAMs positions in pamDB are these:
                        #     NGG  CCN TTTN  NAAA
                        #     ^      ^    ^  ^
                        if motif.extends5 # we branch here because we want to iterate chunk closer to PAM first
                            
                            # PAM10110
                            #   ^ ^  ^ 
                            # PAM10101
                            # PAM11010
                            rve_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[1]])
                            rve_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[2]])
                            sort!(rve_pos_r)
                            err_drift = ch[2] - ch[1] - 1
                            for p in rve_pos_l
                                # find closest chunk to the right
                                stop_l = p + sum(g_chunks_len_rev[ch[1]:(ch[2]-1)]) - 1 - err_drift
                                idx = searchsortedfirst(rve_pos_r, stop_l)
                                if (idx <= length(rve_pos_r) && 
                                    ((rve_pos_r[idx] - stop_l) <= (err_drift * 2 + 1)))

                                    # now find the left side PAM
                                    err_pam = ch[1] - 1
                                    start_nopam = p - sum(g_chunks_len_rev[1:(ch[1]-1)]) + 1 + pam_len + err_pam
                                    # find last value <= to p
                                    idx = searchsortedlast(pam_rve, start_nopam)
                                    # PAM+20-211010
                                    #   ^       ^
                                    while (idx > 0 && length(pam_rve) <= idx && 
                                        (start_nopam - pam_rve[idx] <= (err_pam * 2 + 1)))

                                        pam00!(guides[i], seq, chrom, pam_rve[idx] - pam_len + 1, true, true, motif, rve_filter, detail_file)
                                        idx -= 1
                                    end
                                end
                            end
                        else
                            # 01101PAM
                            # ^  ^ ^ 
                            # 10101PAM
                            # 01011PAM
                            rve_pos_l = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[1]])
                            sort!(rve_pos_l)
                            rve_pos_r = mapreduce(x -> locateall(x, fmi), vcat, g_chunks_rev[ch[2]])
                            err_drift = ch[2] - ch[1] - 1
                            for p in rve_pos_r
                                # find closest chunk to the left
                                start_r = p - sum(g_chunks_len_rev[(ch[1] + 1):(ch[2] - 1)]) + 1 + err_drift
                                idx = searchsortedlast(rve_pos_l, start_r)
                                if idx > 0 && length(rve_pos_l) <= idx && ((start_r - rve_pos_l[idx]) <= (err_drift * 2 + 1))
                                    # now we find closes upstream PAM
                                    err_pam = length(g_chunks_len_rev) - ch[2]
                                    stop_r = p + sum(g_chunks_len_rev[ch[2]:end]) - pam_len - 1 - err_pam
                                    # 01011-20+2PAM
                                    #     ^     ^
                                    idx = searchsortedfirst(pam_rve, stop_r)
                                    while idx <= length(pam_rve) && (pam_rve[idx] - stop_r <= (err_pam * 2 + 1))
                                        pam00!(guides[i], seq, chrom, pam_rve[idx] + pam_len - 1, false, true, motif, fwd_filter, detail_file)
                                        idx += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end

        end

        close(detail_file)
    end
    close(ref)
    cleanup_detail(output_file)
    return
end