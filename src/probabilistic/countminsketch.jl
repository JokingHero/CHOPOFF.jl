"""
    CMSketch{T<:Unsigned}(length, ntables)

Constructs a count-min sketch for memory-friendly counting of hashable objects.
The count returned by a sketch is sometimes higher than the true count, never lower.
The CMSketch does not overflow, - each cell is maximally `typemax(T)`.

# Arguments
* `len`: The number of elements in each hash table
* `ntables`: Number of tables

# Examples
```julia-repl
julia> sketch = CMSketch(5000, 4);
julia> push!(sketch, "hello"); # increment by one
julia> add!(sketch, "hello", eltype(sketch)(5)); # increment by 5
julia> sketch["hello"]
6
```
"""
struct CMSketch{T<:Unsigned}
    len::Int # Cached for speed
    width::Int # Cached for speed
    matrix::Matrix{T}

    function CMSketch{T}(len, ntables) where {T<:Unsigned}
        if len < 1 || ntables < 2
            throw(ArgumentError("Must have len ≥ 1 ntables ≥ 2"))
        end
        return new(len, ntables, zeros(T, (Int(len), Int(ntables))))
    end
end

CMSketch(len, ntables) = CMSketch{UInt8}(len, ntables)

function Base.:(==)(x::CMSketch{T}, y::CMSketch{T}) where {T}
    if x.len != y.len || x.width != y.width
        return false
    end
    return all(i == j for (i,j) in zip(x.matrix, y.matrix))
end

Base.:(==)(x::CMSketch{T1}, y::CMSketch{T2}) where {T1, T2} = false

function Base.show(io::IO, sketch::CMSketch{T}) where {T}
    print(io, "CMSketch{$T}", size(sketch.matrix))
end

index(len, h) = reinterpret(Int, Core.Intrinsics.urem_int(h, reinterpret(UInt64, len))) + 1

@inline function increment!(sketch::CMSketch{T}, h::UInt64, table::Int, count::T) where {T}
    @inbounds existing = sketch.matrix[index(sketch.len, h), table]
    @inbounds sketch.matrix[index(sketch.len, h), table] = safeadd(existing, count)
    return nothing
end

"""
    add!(sketch::CMSketch, val, count)

Add `count` number of `val` to the sketch. For increased speed, let `count` be
of the same type as `eltype(sketch)`.

# Examples
```julia-repl
julia> sketch = CMSketch(1 << 24, 4);
julia> add!(sketch, "hello", eltype(sketch)(5));
julia> sketch["hello"]
5
```
"""
function add!(sketch::CMSketch, x, count)
    # Do not allow negative additions or a count higher than typemax(T)
    # This will screw up saturating arithmetic and guaranteed lower bound.
    count = convert(eltype(sketch), count)
    h = hash(x) # initial hash if it's expensive
    increment!(sketch, h, 1, count)
    for ntable in 2:sketch.width
        #base = (0x93774c8a392b33bb % UInt) * 20
        h = hash(h, reinterpret(UInt64, ntable)) #⊻ base
        increment!(sketch, h, ntable, count)
    end
    return sketch
end

"""
    push!(sketch::CMSketch, val)

Add `val` to the sketch once.
# Examples
```julia-repl
julia> sketch = CMSketch(1 << 24, 4);
julia> push!(sketch, "hello");
julia> sketch["hello"]
1
```
"""
Base.push!(sketch::CMSketch, x) = add!(sketch, x, one(eltype(sketch)))

function Base.append!(sketch::CMSketch, args...)
    for i in vcat(args...)
        push!(sketch, i)
    end
end

function add_guides!(sketch::CMSketch, guides::Vector{LongDNASeq})
    for g in guides
        push!(sketch, unsigned(DNAMer(g)))
    end
end

"""
    haskey(sketch::CMSketch)

Check if sketch[val] > 0.
"""
Base.haskey(sketch::CMSketch, x) = sketch[x] > 0
Base.eltype(sketch::CMSketch{T}) where {T} = T
Base.size(sketch::CMSketch) = (sketch.len, sketch.width)
Base.sizeof(sketch::CMSketch) = 16 + sizeof(sketch.matrix)

"""
    empty!(sketch::CMSketch)

Reset counts of all items to zero, returning the sketch to initial state.
"""
function Base.empty!(sketch::CMSketch)
    fill!(sketch.matrix, zero(eltype(sketch)))
    return sketch
end

"""
    isempty(sketch::CMSketch)

Check if no items have been added to the sketch.
"""
Base.isempty(sketch::CMSketch) = all(i == zero(eltype(sketch)) for i in sketch.matrix[:,1])

function Base.copy!(dst::CMSketch{T}, src::CMSketch{T}) where {T}
    if dst.len != src.len || dst.width != src.width
        throw(ArgumentError("Sketches must have same dimensions"))
    end
    unsafe_copyto!(dst.matrix, 1, src.matrix, 1, src.len * src.width)
    return dst
end

function Base.copy(sketch::CMSketch)
    newsketch = typeof(sketch)(sketch.len, sketch.width)
    return copy!(newsketch, sketch)
end

"""
    +(x::CMSketch, y::CMSketch)

Add two count-min sketches together. Will not work if `x` and `y` do not share
parameters `T`, `length` and `width`. The result will be a sketch with the summed
counts of the two input sketches.

# Examples
```
julia> x, y = CMSketch(1000, 4), CMSketch(1000, 4);

julia> add!(x, "hello", 4); add!(y, "hello", 19);

julia> z = x + y; Int(z["hello"])
23
```
"""
function Base.:+(x::CMSketch{T}, y::CMSketch{T}) where {T}
    if x.len != y.len || x.width != y.width
        throw(ArgumentError("Sketches must have same dimensions"))
    end
    summed = copy(x)
    for i in 1:(x.len * x.width)
        @inbounds summed.matrix[i] = safeadd(summed.matrix[i], y.matrix[i])
    end
    return summed
end


"
Calculate at which percent database is filled, preferably 
should be below 0.8.

Estimate the probability of miscounting an element in the sketch.
"
function fillrate(sketch::CMSketch)
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
function zerorate(sketch::CMSketch)
    return sum(sketch.matrix .> zero(eltype(sketch))) / (sketch.len * sketch.width)
end

"""
    getindex(sketch::CMSketch, item)

Get the estimated count of `item`. This is never underestimated, but may be
overestimated.
"""
function Base.getindex(sketch::CMSketch, x)
    h = hash(x) # initial hash if it's expensive
    #base = (0x93774c8a392b33bb % UInt) * 20
    @inbounds count = sketch.matrix[index(sketch.len, h), 1]
    for ntable in 2:sketch.width
        h = hash(h, reinterpret(UInt64, ntable)) #⊻ base
        @inbounds m = sketch.matrix[index(sketch.len, h), ntable]
        count = min(count, m)
    end
    return count
end