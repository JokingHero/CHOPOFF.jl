struct SketchDB{T<:Unsigned}
    sketch::CountMinSketch{T}
    kmers::IdDict{String, Int}
    dbi::DBInfo
end


"
Calculate at which percent database is filled, preferably 
should be below 0.8.
"
function fillrate(sketch::CountMinSketch)
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


function makeemptysketch(
    est_len::Int,
    max_count::Int,
    probability_of_error = 0.01)

    max_count_type = smallestutype(unsigned(max_count))
    # we compute how many hashing functions we need
    depth = ceil(log(1/probability_of_error))
    # estimate our error E based on the max_len
    E = 1/(est_len)
    width = ceil(Base.ℯ/E)
    sketch = CountMinSketch{max_count_type}(width, depth)
    return sketch
end


"
Build Count-Min-Sketch and kmer count database.
"
function buildsketchDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    probability_of_error::Float64 = 0.01; 
    max_count::Int = 255)

    dbi = DBInfo(genomepath, name, motif)

    # first we measure how many unique guides there are
    @info "Unique guide count estimations..."
    hll = HyperLogLog{18}()
    gatherofftargets!(hll, dbi)

    # next we adjust the estimated length +5% for error
    est_len = Int(ceil(1.05*length(hll)))
    @info "We estimate $est_len unique guides."
    @info "Constructing sketch database..."
    sketch = makeemptysketch(est_len, max_count, probability_of_error)
    gatherofftargets!(sketch, dbi)

    @info "Constructing kmer sketch..."
    kmers = IdDict{String, Int}()
    gatherofftargets!(kmers, dbi)

    db = SketchDB(sketch, kmers, dbi)
    save(db, joinpath(storagedir, "sketchDB.bin"))
    @info "Finished constructing sketchDB in " * storagedir
    @info "Database fillrate is " * string(round(fillrate(db.sketch); digits = 3))
    return storagedir
end


"
For each of the input guides will check in the sketchDB 
what value Count-Min Sketch has for it in the first, and 
in the second column will output kmer sum. Results 
are estimations of offtarget counts in the genome.

If CMS column is 0, it is guaranteed this guide has 
no 0-distance off-targets in the genome!
"
function searchsketchDB(
    storagedir::String,
    guides::Vector{LongDNASeq})

    sdb = load(joinpath(storagedir, "sketchDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    kmer_len = length(first(keys(sdb.kmers)))
    # TODO check that seq is in line with motif

    res = zeros(Float64, length(guides_), 2)
    for (i, s) in enumerate(guides_)
        res[i, 1] += sdb.sketch[string(s)] # 0 distance

        counts = Vector()
        for j in 1:(length(guides_[i]) - kmer_len) # 0-max_dist
            push!(counts, get(sdb.kmers, string(guides_[i][j:(j + kmer_len - 1)]), 0))
        end
        res[i, 2] = maximum(counts) # sum is worst!, median best
    end

    res = DataFrame(res)
    col_d = ["CMS", "Kmer"]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end