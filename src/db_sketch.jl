struct SketchDB{T<:Unsigned}
    sketch::CountMinSketch{T}
    dbi::DBInfo
end


"
Calculate at which percent database is filled, preferably 
should be below 0.8.

Estimate the probability of miscounting an element in the sketch.
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


"
What % of the database contains zeros.
"
function zerorate(sketch::CountMinSketch)
    return sum(sketch.matrix .> zero(eltype(sketch))) / (sketch.len * sketch.width)
end


function makeemptysketch(
    est_len::Int,
    max_count::Int,
    probability_of_error = 0.001)

    max_count_type = smallestutype(unsigned(max_count))
    # we compute how many hashing functions we need
    depth = ceil(log(1/probability_of_error))
    # estimate our error E based on the max_len
    width = ceil(Base.ℯ / (1/(est_len)))
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
    probability_of_error::Float64 = 0.001; 
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
    
    db = SketchDB(sketch, dbi)
    save(db, joinpath(storagedir, "sketchDB.bin"))
    @info "Finished constructing sketchDB in " * storagedir
    @info "Estimated probability of miscounting an element in the sketch is " * string(round(fillrate(db.sketch); digits = 6))
    @info "Database empty field fraction is " * string(round(zerorate(db.sketch); digits = 3))
    @info "Database size is:\n width -> " * string(db.sketch.width) *
        "\n length -> " * string(db.sketch.len) *
        "\n consuming: " * string(round((sizeof(db) * 1e-6); digits = 3)) * " mb of space."
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
    guides::Vector{LongDNASeq},
    dist::Int = 2)

    sdb = load(joinpath(storagedir, "sketchDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), dist + 1)
    for (i, s) in enumerate(guides_)
        s = string(s)
        res[i, 1] += sdb.sketch[s] # 0 distance
        for d in 1:dist
            res[i, d + 1] = sum([sdb.sketch[sd] for sd in comb_of_d(s, d)])
        end
    end

    res = DataFrame(res)
    col_d = [Symbol("D$i") for i in 0:dist]
    rename!(res, col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end