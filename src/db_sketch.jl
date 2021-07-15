struct SketchDB{T<:Unsigned}
    sketch::CountMinSketch{T}
    dbi::DBInfo
end


"
What % of the database contains zeros.
"
function zerorate(sketch::CountMinSketch)
    return sum(sketch.matrix .> zero(eltype(sketch))) / (sketch.len * sketch.width)
end


function add_guides!(sketch::HyperLogLog, guides::Vector{LongDNASeq})
    for g in guides
        if isambiguous(g)
            g = expand_ambiguous(g)
            for gi in g
                push!(sketch, gi)
            end
        else
            push!(sketch, g)
        end
    end
end


function add_guides!(sketch::CountMinSketch, guides::Vector{LongDNASeq})
    for g in guides
        push!(sketch, g)
    end
end


function makeemptysketch(
    est_len::Int,
    max_count::Int,
    probability_of_error = 0.001,
    error_size = 3)

    max_count_type = smallestutype(unsigned(max_count))
    # we compute how many hashing functions we need
    depth = max(ceil(log(1/probability_of_error)), 2)
    # estimate our error E based on the max_len
    width = max(ceil(Base.ℯ / (error_size/(est_len))), 2)
    sketch = CountMinSketch{max_count_type}(width, depth)
    return sketch
end


"""
`build_sketchDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    probability_of_error::Float64 = 0.001; 
    max_count::Int = 255)`

Build Count-Min-Sketch database that can never under-estimate the counts of
off-targets.

# Arguments

`name` - Your name for this database.  
`genomepath` - Path to the genome in either .2bit or .fa format.
`motif`- See Motif type.  
`storagedir`- Where should be the output folder.  
`probability_of_error`- Transforms into number of used hash functions. 
    This is error level for each individual off-target search. 
    This error will be propagated when doing searches for larger distances and therefore 
    will not be perfectly reflected in the results.
`error_size` - Transforms into length of the table, it is epsilon indicating what is desired
    error of the estimated count we can afford.

`max_count` - What is the maximum count for each specific off-target in the genome. We
    normally care only for guides that have low off-target count, therefore we can keep this
    value also low, this saves space and memory usage.
"""
function build_sketchDB(
    name::String, 
    genomepath::String, 
    motif::Motif,
    storagedir::String,
    probability_of_error::Float64 = 0.001,
    error_size::Int = 3;
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
    sketch = makeemptysketch(est_len, max_count, probability_of_error, error_size)
    gatherofftargets!(sketch, dbi)
    
    db = SketchDB(sketch, dbi)
    save(db, joinpath(storagedir, "sketchDB.bin"))
    @info "Finished constructing sketchDB in " * storagedir
    @info "Estimated probability of miscounting an element in the sketch is " * string(round(fprof(db.sketch); digits = 6))
    @info "Database empty field fraction is " * string(round(zerorate(db.sketch); digits = 3))
    @info "Database size is:\n width -> " * string(db.sketch.width) *
        "\n length -> " * string(db.sketch.len) *
        "\n consuming: " * string(round((filesize(joinpath(storagedir, "sketchDB.bin")) * 1e-6); digits = 3)) * 
        " mb of disk space."
    return storagedir
end


"""
`search_sketchDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)`

For each of the input guides will check in the sketchDB 
what value Count-Min Sketch has for it in the first, and 
in the second column will output kmer sum. Results 
are estimations of offtarget counts in the genome.

If CMS column is 0, it is guaranteed this guide has 
no 0-distance off-targets in the genome!

Also, maximum count for each off-target in the database is capped,
to `max_count` specified during building of sketchDB. This means 
that counts larger than `max_count` (default 255) are no longer
estimating correctly, but we know at least 255 off-targets exist in 
the genome. Its likely you would not care for those guides anyhow.

In principle, you can use `distance` of any size, but for each increasing distance value
we have recurent solution that slows down significantly from distance 3 and upwards. These estimations
should be used normally just as prefiltering step for the guideRNA's, therefore checking just distance 0,
and distance 1 off-targets should be enough to rank the guides. For increasing distances partial errors in
estimations are stacking and overwhelm real estimations quickly.

## Result

Returns .csv table with the following columns:

D0, D1 ... D`distance` - these are estimated counts of off-targets for given guideRNA, 
each column for increasing `distance`, these is also a sum of the corresponding DN and DB columns

DN1, DN2, ... DN`distance`- these values are estimated counts of off-targets for given guideRNA,
each column for increasing `distance`, but these values are for guaranteed "within" distance off-targets

DB1, DB2, ... DB`distance` - same as above, but these values are representing "borderline cases" where 
genomic extension is not known and only asummed the worse case. These values are assuming worst case, 
and therefore are overestimating extensively with increasing `distance`.
"""
function search_sketchDB(
    storagedir::String,
    guides::Vector{LongDNASeq},
    dist::Int = 1)

    if any(isambiguous.(guides))
        throw("Ambiguous bases are not allowed in guide queries.")
    end

    sdb = load(joinpath(storagedir, "sketchDB.bin"))
    guides_ = copy(guides)
    # reverse guides so that PAM is always on the left
    if sdb.dbi.motif.extends5
        guides_ = reverse.(guides_)
    end

    # TODO check that seq is in line with motif
    res = zeros(Int, length(guides_), (dist + 1) * 3 - 2)
    for (i, s) in enumerate(guides_)
        res[i, 1] += sdb.sketch[s]
        for d in 1:dist
            norm_d, border_d = comb_of_d(string(s), d)
            norm_d_res = ThreadsX.sum(sdb.sketch[LongDNASeq(sd)] for sd in norm_d)
            border_d_res = ThreadsX.sum(sdb.sketch[LongDNASeq(sd)] for sd in border_d)
            res[i, d + 1] = norm_d_res + border_d_res
            res[i, dist + d + 1] = norm_d_res
            res[i, dist * 2 + d + 1] = border_d_res
        end
    end

    res = DataFrame(res, :auto)
    col_d = [Symbol("D$i") for i in 0:dist]
    all_col_d = vcat(col_d, [Symbol("DN$i") for i in 1:dist], [Symbol("DB$i") for i in 1:dist])
    rename!(res, all_col_d)
    res.guide = guides
    sort!(res, col_d)
    return res
end