"
Find all instances of pat inside seq, disregard seq ambiguous bases.
Restrict seq to subset of start:stop.
"
function Base.findall(pat::BioSequence, seq::BioSequence,
    start::Integer = 1, stop::Integer = lastindex(seq))

    res = Vector{UnitRange{Int64}}()
    m = length(pat)
    n = length(seq)
    stop_ = min(stop, n) - m
    s::Int = max(start - 1, 0)

    if m == 0  # empty query
        return nothing
    end

    while s ≤ stop_
        if isambiguous(seq[s+m])
            s += m
        elseif iscompatible(pat[m], seq[s+m])
            i = m - 1
            while i > 0
                if isambiguous(seq[s+i]) || !iscompatible(pat[i], seq[s+i])
                    break
                end
                i -= 1
            end
            if i == 0
                push!(res, (s+1):(s+length(pat))) #found
            end
            s += 1
        else
            s += 1
        end
    end

    return res  # not found
end

"
Removes PAM from the seq.
"
function removepam(seq::LongDNASeq, pam::Vector{UnitRange{<:Integer}})
    if length(pam) == 1
        seq = copy(seq)
        deleteat!(seq, pam[1])
    else
        pam = vcat([collect(x) for x in pam]...)
        seq = LongDNASeq(string(seq)[pam])
    end
    return seq
end


function pushguides!(
    chrom::K,
    query::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    max_dist::Int,
    reverse_comp::Bool,
    output::T) where {T<:Union{CountMinSketch, HyperLogLog, Vector{String}}, K<:BioSequence}

    if length(query) != 0
        for x in findall(query, chrom)
            guide = chrom[x]
            guide = removepam(guide, pam_loci)
            if reverse_comp
                reverse_complement!(guide)
            end

            push!(output, string(guide)) # 0 distance
        end
    end
    return output
end


function pushguides!(
    chrom::K,
    query::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    max_dist::Int,
    reverse_comp::Bool,
    output::Dict{String, Int}) where {K<:BioSequence}

    g_len = length(query) - length(collect(vcat(pam_loci...)))
    kmer_size = Int(floor(g_len / (max_dist + 1)))

    if length(query) != 0
        for x in findall(query, chrom)
            guide = chrom[x]
            guide = removepam(guide, pam_loci)
            if reverse_comp
                reverse_complement!(guide)
            end

            # kmer treatment - just naive count
            for i in 1:g_len-kmer_size
                kmer = string(guide[i:(i + kmer_size - 1)])
                output[kmer] = get(output, kmer, 0) + 1
            end
        end
    end
    return output
end


function findofftargets!(
    genome::String,
    motif::Motif,
    max_dist::Int,
    output::T) where T<:Union{CountMinSketch, HyperLogLog, Vector{String}, Dict}

    ext = extension(genome)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        throw(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(genome, "r")
    reader = is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)

    for record in reader
        chrom = is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)

        if FASTA.hasidentifier(record) # TODO 2bit?!
            println("Working on: ", FASTA.identifier(record))
        end

        pushguides!(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, output)
        pushguides!(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, output)
    end

    close(ref)
    return
end


function makeemptysketch(
    est_len::Int,
    max_count::Unsigned,
    probability_of_error = 0.01)

    max_count_type = smallestutype(max_count)
    # we compute how many hashing functions we need
    depth = ceil(log(1/probability_of_error))
    # estimate our error E based on the max_len
    # we add 0.5% error of estimation for HLL
    E = 1/(est_len)
    width = ceil(Base.ℯ/E)
    sketch = CountMinSketch{max_count_type}(width, depth)
    return sketch
end


function findofftargets(
    genome::String,
    motif::Motif,
    max_dist::Int,
    max_count = 255)

    # first we measure how many unique guides there are
    println("Unique guide count estimations...")
    hll = HyperLogLog{18}()
    findofftargets!(genome, motif, max_dist, hll)

    # next we adjust the estimated length +5% for error
    est_len = Int(ceil(5.05*length(hll)))
    println("HLL estimations are: ", est_len)

    println("Constructing sketch database...")
    sketch = makeemptysketch(est_len, unsigned(max_count))
    findofftargets!(genome, motif, max_dist, sketch)

    println("Constructing kmer sketch...")
    kmers = Dict{String, Int}() # TODO IdDict would be better
    findofftargets!(genome, motif, max_dist, kmers)

    println("Done!")
    db = SketchDB(sketch, kmers, motif, max_dist)
    return db
end


function estimate(
    db::SketchDB,
    seq::Vector{LongDNASeq}, # asume seq has no PAM and is 5'-3' direction
    max_dist = db.max_dist)

    max_dist = max_dist > db.max_dist ? db.max_dist : max_dist
    kmer_len = length(first(keys(db.kmers)))
    # TODO check that seq is in line with motif

    res = zeros(Int, length(seq), 2)
    for (i, s) in enumerate(seq)
        res[i, 1] += db.sketch[string(s)] # 0 distance

        counts = Vector()
        for j in 1:(length(seq[i]) - kmer_len) # 0-max_dist
            push!(counts, get(db.kmers, string(seq[i][j:(j + kmer_len - 1)]), 0))
        end
        res[i, 2] = mean(counts)
    end

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:1]
    rename!(res, col_d)
    res.guide = seq
    sort!(res, col_d)
    return res
end


function fprof_fixed(sketch::CountMinSketch)
    rate = 1
    for col in 1:sketch.width
        full_in_row = 0
        for row in 1:sketch.len
            full_in_row += sketch.matrix[row, col] > zero(eltype(sketch))
        end
        rate *= full_in_row / sketch.len
    end
    return rate
end


function fillrate(db::SketchDB)
    return fprof_fixed(db.sketch)
end



function countofftargets(
    chrom::K,
    query::LongDNASeq,
    pam_loci::Vector{UnitRange{<:Integer}},
    max_dist::Int,
    reverse_comp::Bool,
    guides::Vector{LongDNASeq}) where {K<:BioSequence}

    output = zeros(Int, length(guides), max_dist + 1)

    if length(query) != 0
        for x in findall(query, chrom)
            guide_ref = chrom[x]
            guide_ref = removepam(guide_ref, pam_loci)
            if reverse_comp
                reverse_complement!(guide_ref)
            end

            for (i, g) in enumerate(guides)
                dist = levenshtein(reverse(g), reverse(guide_ref), max_dist)
                if dist <= max_dist
                    output[i, dist + 1] += 1
                end
            end
        end
    end
    return output
end


function findofftargets(
    genome::String,
    motif::Motif, # has to be motive including the max_dist extensions on the 3'
    max_dist::Int,
    g::Array{LongDNASeq, 1}) # vector of guides in 1-20 (NGG) config)

    ext = extension(genome)
    if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        is_fa = true
    elseif (ext == ".2bit")
        is_fa = false
    else
        throw(
        "Wrong extension of the genome.",
        "We can parse only .fa/.fna/.fasta or .2bit references.")
    end

    ref = open(genome, "r")
    reader = is_fa ? FASTA.Reader(ref) : TwoBit.Reader(ref)
    reader = collect(enumerate(reader))
    shuffle!(reader)

    res = zeros(Int, length(g), max_dist + 1)
    res_vec = [res for r in 1:length(reader)]
    Threads.mapi()
    for (i, record) in reader #TODO use  instead
        chrom = is_fa ? FASTA.sequence(record) : TwoBit.sequence(record)

        if FASTA.hasidentifier(record) # TODO 2bit?!
            println("Working on: ", FASTA.identifier(record))
        end

        res_vec[i] += countofftargets(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, g)
        res_vec[i] += countofftargets(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, g)
    end
    close(ref)

    res = reduce(+, res_vec)
    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:max_dist]
    rename!(res, col_d)
    res.guide = [string(gi) for gi in g]
    sort!(res, col_d)
    return res

end
