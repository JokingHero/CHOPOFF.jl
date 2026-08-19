"""
```
Motif(
    alias::String, 
    fwdmotif::String, 
    fwdpam::String, 
    forward_strand::Bool = true, 
    reverse_strand::Bool = true, 
    distance::Int = 4, 
    extends5::Bool = true,
    ambig_max::Int = 5)
```

`Motif(alias::String)`


Motif defines what we search on the genome,
what can be identified as an off-target.


# Arguments
`alias` - alias of the motif for easier identification e.g. Cas9

`fwdmotif` - Motif that indicates where is PAM inside `fwdpam`.
    For example for Cas9 it is 20*N + XXX: NNNNNNNNNNNNNNNNNNNNXXX

`fwdpam`   - Motif in 5'-3' that will be matched on the reference (without the X).
             For example for Cas9 it is 20*X + NGG:
             XXXXXXXXXXXXXXXXXXXXNGG

`forward`  - If false will not match to the forward reference strand.

`reverse`  - If false will not match to the reverse reference strand.

`distance` - How many extra nucleotides are needed for a search? This
             will indicate within what distance we can search for off-targets.
             When we don't have those bases we use DNA_Gap.
             
`extend5`  - Defines how off-targets will be aligned to the guides and where
             extra nucleotides will be added for alignment within distance. Whether
             to extend in the 5' and 3' direction. Cas9 is extend5 = true.

`ambig_max`- How many ambiguous bases are allowed in the pattern?

```
Example for Cas9 where we want to search for off-targets within distance of 4:
  alias:    Cas9
  fwdmotif: NNNNNNNNNNNNNNNNNNNNXXX
  fwdpam:   XXXXXXXXXXXXXXXXXXXXNGG
  forward:  true
  reverse:  true
  distance: 4
  extend5:  true
  ambig_max:5 
```

Alignments will be performed from opposite to the extension direction (which is defined by extend5).

# Examples
```jldoctest
julia> Motif("Cas9")
Alias: Cas9
Maximum search distance: 4
Number of allowed ambigous bp: 0
20N-NGG

julia> Motif("test name", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 4, true, 5)
Alias: test name
Maximum search distance: 4
Number of allowed ambigous bp: 5
20N-NGG
```
"""
struct Motif
    alias::String
    fwd::LongDNA{4}
    rve::LongDNA{4}
    pam_loci_fwd::UnitRange{<:Integer}
    pam_loci_rve::UnitRange{<:Integer}
    distance::Int
    extends5::Bool
    ambig_max::Int
end

Base.:(==)(x::Motif, y::Motif) = x.fwd == y.fwd && 
    x.rve == y.rve &&
    x.pam_loci_fwd == y.pam_loci_fwd && 
    x.pam_loci_rve == y.pam_loci_rve &&
    x.distance == y.distance &&
    x.extends5 == y.extends5 &&
    x.ambig_max == y.ambig_max


function notX(s1, s2, x = 'X')
    if s1 != x
        return s1
    end

    if s2 != x
        return s2
    else
        throw("Both letters are X.")
    end
end


function combinestrings(s1::String, s2::String, rule::Function = notX)
    if (length(s1) != length(s2))
        error("Unequal lengths.")
    end
    return join([rule(s1[i], s2[i]) for i in eachindex(s1)])
end


"
Removes PAM from the seq.
"
function removepam(seq::LongDNA{4}, pam::UnitRange{<:Integer})
    x = copy(seq)
    deleteat!(x, pam)
    return x
end


function Motif(alias::String,
    fwdmotif::String, fwdpam::String,
    forward_strand::Bool = true, reverse_strand::Bool = true,
    distance::Int = 4, extends5::Bool = true, ambig_max::Int = 0)
    if length(fwdmotif) != length(fwdpam)
        throw("fwd_motif and fwd_pam have to have the same length!")
    end
    merge = combinestrings(fwdmotif, fwdpam)
    pam_matches = findall(r"[^X]+", fwdpam)
    length(pam_matches) <= 1 ||
        error("PAM bases must form one contiguous block.")
    pam_loci = isempty(pam_matches) ? (1:0) : only(pam_matches)

    if forward_strand
        pam_loci_fwd = pam_loci
        fwd = LongDNA{4}(merge)
    else
        pam_loci_fwd = 1:0
        fwd = LongDNA{4}("")
    end

    if reverse_strand
        pam_loci_rve = isempty(pam_loci) ? (1:0) : (
            (length(fwdpam) - last(pam_loci) + 1):
            (length(fwdpam) - first(pam_loci) + 1))
        rve = reverse_complement(LongDNA{4}(merge))
    else
        pam_loci_rve = 1:0
        rve = LongDNA{4}("")
    end

    return Motif(alias, fwd, rve, pam_loci_fwd, pam_loci_rve, distance, extends5, ambig_max)
end


"""
`length_noPAM(motif::Motif)`

Calculate what is the length of the motif, without extension, and without PAM.
Effectively, size of the gRNA.

# Examples
```jldoctest
julia> length_noPAM(Motif("Cas9"))
20
```
"""
function length_noPAM(motif::Motif)
    fwd_len = length(motif.fwd) - length(motif.pam_loci_fwd)
    rve_len = length(motif.rve) - length(motif.pam_loci_rve)
    return max(fwd_len, rve_len)
end


import Base.length 
"""
`length(motif::Motif)`

Length of the motif with PAM, without extension.

# Examples
```jldoctest
julia> length(Motif("Cas9"))
23
```
"""
function length(motif::Motif)
    return max(length(motif.fwd), length(motif.rve))
end


"""
`setambig(motif::Motif, ambig::Int)`

Set the ambiguity limit across the complete guide and PAM window. Ambiguous
bases in the edit-distance extension are not counted.

# Examples
```jldoctest
julia> setambig(Motif("Cas9"), 15)
Alias: Cas9
Maximum search distance: 4
Number of allowed ambigous bp: 15
20N-NGG
```
"""
function setambig(motif::Motif, ambig::Int)
    return Motif(motif.alias, motif.fwd, motif.rve, 
        motif.pam_loci_fwd, motif.pam_loci_rve, 
        motif.distance, motif.extends5, ambig)
end


"""
`setdist(motif::Motif, distance::Int)`

Set the distance (maximum value of allowed mismatches, deletion, insertions) 
that are allowed during alignment.

# Examples
```jldoctest
julia> setdist(Motif("Cas9"), 2)
Alias: Cas9
Maximum search distance: 2
Number of allowed ambigous bp: 0
20N-NGG
```
"""
function setdist(motif::Motif, distance::Int)
    return Motif(motif.alias, motif.fwd, motif.rve, 
        motif.pam_loci_fwd, motif.pam_loci_rve, 
        distance, motif.extends5, motif.ambig_max)
end


const motif_db = Dict(
    "test" => Motif("test",
                    "NNNX",
                    "XXXG", true, true, 2, true, 0),
    "Cas9" => Motif("Cas9",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 4, true, 0),
    "Cas9_NNN" => Motif("Cas9_NNN",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNNN", true, true, 4, true, 0),
    "Cas9_NGA" => Motif("Cas9_NGA",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNGA", true, true, 4, true, 0),
    "Cas9_NAG" => Motif("Cas9_NAG",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNAG", true, true, 4, true, 0),
    "Cas9_NNGT" => Motif("Cas9_NNGT",
                    "NNNNNNNNNNNNNNNNNNNNXXXX",
                    "XXXXXXXXXXXXXXXXXXXXNNGT", true, true, 4, true, 0),
    "Cas9_NNAGAA" => Motif("Cas9_NNAGAA",
                    "NNNNNNNNNNNNNNNNNNNNXXXXXX",
                    "XXXXXXXXXXXXXXXXXXXXNNAGAA", true, true, 4, true, 0),
    "Cas9_NGGNG" => Motif("Cas9_NGGNG",
                    "NNNNNNNNNNNNNNNNNNNNXXXXX",
                    "XXXXXXXXXXXXXXXXXXXXNGGNG", true, true, 4, true, 0),
    "Cas12a" => Motif("Cas12a",
                    "XXXXNNNNNNNNNNNNNNNNNNNNN",
                    "TTTVXXXXXXXXXXXXXXXXXXXXX", true, true, 4, false, 0),
    "CasX" => Motif("CasX",
                    "XXXXNNNNNNNNNNNNNNNNNNNN",
                    "TTCNXXXXXXXXXXXXXXXXXXXX", true, true, 4, false, 0),
    "hfCas12Max" => Motif("hfCas12Max",
                    "XXXNNNNNNNNNNNNNNNNNNNNNNN",
                    "TNNXXXXXXXXXXXXXXXXXXXXXXX", true, true, 4, false, 0)
    )


function Motif(alias::String; distance::Int = 4, ambig_max::Int = 0)
    motif = motif_db[alias]
    return setambig(setdist(motif, distance), ambig_max)
end


# all sequences will be ready to search on forward strand here
# offtarget should be reversed if we are extend5 true
# COV_EXCL_START
function appendPAM_forward(offtarget::LongDNA{4}, motif::Motif)
    if motif.extends5
        return offtarget * motif.fwd[motif.pam_loci_fwd]
    end
    return motif.fwd[motif.pam_loci_fwd] * offtarget
end
# COV_EXCL_STOP


function visualize_motif(motif::Motif; separator::String = "_")
    reverse_only = isempty(motif.fwd)
    pattern, pam_loci = reverse_only ?
        (motif.rve, motif.pam_loci_rve) :
        (motif.fwd, motif.pam_loci_fwd)
    isempty(pam_loci) && return string(length_noPAM(motif)) * "N"
    pam = pattern[pam_loci]
    reverse_only && (pam = reverse_complement(pam))
    s = string(pam)
    if motif.extends5
        s = string(length_noPAM(motif)) * "N" * separator * s
    else
        s = s * separator * string(length_noPAM(motif)) * "N"
    end
    return s
end


function display_motif(motif::Motif)
    return "Alias: " * motif.alias * "\n" * 
        "Maximum search distance: " * string(motif.distance) * "\n" * 
        "Number of allowed ambigous bp: " * string(motif.ambig_max) * "\n" * 
        visualize_motif(motif; separator = "-")
end

Base.show(io::IO, ::MIME"text/plain", f::Motif) = print(io, display_motif(f))
Base.print(io::IO, f::Motif) = print(io, f.alias * "_" * visualize_motif(f) * "_maxDist" * string(f.distance) * "_maxAmbig" * string(f.ambig_max))
