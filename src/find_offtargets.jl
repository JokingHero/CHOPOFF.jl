
"
Removes PAM from the seq.
"
function removepam(seq::LongDNASeq, pam::Vector{UnitRange{<:Integer}})
    if length(pam) == 1
        seq = copy(seq)
        deleteat!(seq, pam[1])
    else
        pam = vcat([collect(x) for x in pam]...)
        seq = string(seq)[pam]
        seq = LongDNASeq(seq)
    end
    return seq
end


function pushdeletes!(
    chrom::K,
    query::ExactSearchQuery{LongSequence{DNAAlphabet{4}}},
    pam_loci::Vector{UnitRange{<:Integer}},
    distance::Int,
    reverse_comp::Bool,
    output::Vector{T}) where {T<:Union{CountMinSketch, HyperLogLog}, K<:BioSequence}

    if length(query.seq) != 0
        start = 1
        guide_len = length(removepam(query.seq, pam_loci))
        deletion_perm = deletion_permutations(guide_len, distance, 1)
        deletion_perm = [sort(vcat(x...)) for x in deletion_perm]
        deletion_perm = [setdiff(collect(1:guide_len), x) for x in deletion_perm]

        while start != nothing
            loci = findfirst(query, chrom, start)
            if loci == nothing
                start = nothing
            else
                guide = removepam(chrom[loci], pam_loci)
                if reverse_comp
                    reverse_complement!(guide)
                end

                if !hasambiguity(guide)
                    push!(output[1], string(guide)) # 0 distance
                    # for the rest of the deletes
                    deletes = unique([string(guide)[x] for x in deletion_perm])
                    for del in deletes
                        push!(output[guide_len - length(del) + 1], del)
                    end
                end
                start = loci.start + 1
            end
        end
    end
    return output
end


function findofftargets!(
    genome::String,
    motif::Motif,
    max_dist::Int,
    output::Vector{T}) where T<:Union{CountMinSketch, HyperLogLog}

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

        pushdeletes!(chrom, motif.fwd, motif.pam_loci_fwd, max_dist, false, output)
        pushdeletes!(chrom, motif.rve, motif.pam_loci_rve, max_dist, true, output)
    end

    close(ref)
    return
end


function makeemptysketches(
    est_len::Vector{Int},
    max_count::Unsigned,
    probability_of_error = 0.01)

    max_count_type = smallestutype(max_count)
    sketches = Vector{CountMinSketch{max_count_type}}() # 1st sketch => 0 distance, and so on
    for max_len in est_len
        # we compute how many hashing functions we need
        depth = ceil(log(1/probability_of_error))
        # estimate our error E based on the max_len
        # we add 0.5% error of estimation for HLL
        E = 1/(max_len)
        width = ceil(Base.ℯ/E)
        push!(sketches, CountMinSketch{max_count_type}(width, depth))
    end
    return sketches
end


function findofftargets(
    genome::String,
    motif::Motif,
    max_dist::Int,
    max_count = 255)
    # first we measure how many unique guides there are
    println("HLL estimations...")
    hll = [HyperLogLog{18}() for i in 1:(max_dist + 1)]
    findofftargets!(genome, motif, max_dist, hll)

    # next we adjust the estimated length +5% for error
    est_len = [Int(ceil(length(x) + 0.05*length(x))) for x in hll]
    println("HLL estimations are: ", est_len)
    sketches = makeemptysketches(est_len, unsigned(max_count))
    println("Constructing database...")
    findofftargets!(genome, motif, max_dist, sketches)
    println("Done!")
    db = SketchDB(sketches, motif)
    return db
end


function estimate(
    db::SketchDB,
    seq::Vector{LongDNASeq}, # asume seq has no PAM and is 5'-3' direction
    max_dist = (length(db.sketches) - 1))

    sketches = db.sketches
    max_dist = max_dist > (length(sketches) - 1) ? (length(sketches) - 1) : max_dist
    # TODO check that seq is in line with motif

    # compute deletes
    guide_len = length(seq[1])
    @assert all([length(x) == guide_len for x in seq])
    deletion_perm = deletion_permutations(guide_len, max_dist, 1)
    deletion_perm = [sort(vcat(x...)) for x in deletion_perm]
    deletion_perm = [setdiff(collect(1:guide_len), x) for x in deletion_perm]

    res = zeros(Int, length(seq), max_dist + 1)
    for (i, s) in enumerate(seq)
        res[i, 1] += sketches[1][string(s)] # 0 distance
        deletes = unique([string(s)[x] for x in deletion_perm])
        for del in deletes # other distances
            col = guide_len - length(del) + 1
            res[i, col] += sketches[col][del]
        end
    end

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:max_dist]
    rename!(res, col_d)
    res.guide = seq
    sort!(res, col_d)
    return res
end
