struct SketchDB{T<:Unsigned}
    sketch::CountMinSketch{T}
    kmers::IdDict{String, Int}
    dbi::DBInfo
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

"
Build Count-Min-Sketch and kmer count databse.
"
function buildsketch(name::String, genomepath::String, motif::Motif,
                     probability_of_error = 0.01, max_count = 255)

    # Save basic gneome information
    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    @info "Unique guide count estimations..."
    hll = HyperLogLog{18}()
    gatherofftargets!(hll, dbi)

    # next we adjust the estimated length +5% for error
    est_len = Int(ceil(1.05*length(hll)))
    @info "HLL estimations are: $est_len"

    @info "Constructing sketch database..."
    sketch = makeemptysketch(est_len, unsigned(max_count))
    gatherofftargets!(sketch, dbi)

    println("Constructing kmer sketch...")
    kmers = IdDict{String, Int}()
    gatherofftargets!(kmers, dbi)

    println("Done!")
    db = SketchDB(sketch, kmers, dbi)
    return db
end


# TODO NOT CHECKED YET TO WORK

function estimate(
    db::SketchDB,
    seq::Vector{LongDNASeq}) # asume seq has no PAM and is 5'-3' direction

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
