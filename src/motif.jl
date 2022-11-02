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
Motif("Cas9", NNNNNNNNNNNNNNNNNNNNNGG, CCNNNNNNNNNNNNNNNNNNNNN, 21:23, 1:3, 3, true, 0)

julia> Motif("test name", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 3, true, 5)
Motif("test name", NNNNNNNNNNNNNNNNNNNNNGG, CCNNNNNNNNNNNNNNNNNNNNN, 21:23, 1:3, 3, true, 5)
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
    distance::Int = 3, extends5::Bool = true, ambig_max::Int = 3)
    if length(fwdmotif) != length(fwdpam)
        throw("fwd_motif and fwd_pam have to have the same length!")
    end
    merge = combinestrings(fwdmotif, fwdpam)

    if forward_strand
        # where is PAM located?
        pam_loci_fwd = findall(r"[^X]+", fwdpam)[1]
        fwd = LongDNA{4}(merge)
    else
        pam_loci_fwd = UnitRange{Int64}()
        fwd = LongDNA{4}("")
    end

    if reverse_strand
        pam_loci_rve = findall(r"[^X]+", reverse(fwdpam))[1]
        rve = reverse_complement(LongDNA{4}(merge))
    else
        pam_loci_rve = UnitRange{Int64}()
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

Set the ambiguity (how many ambiguous bases are allowed, not counting PAM, not counting extension) level for `motif`.

# Examples
```jldoctest
julia> setambig(Motif("Cas9"), 15)
Motif("Cas9", NNNNNNNNNNNNNNNNNNNNNGG, CCNNNNNNNNNNNNNNNNNNNNN, 21:23, 1:3, 3, true, 15)
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
julia> setdist(Motif("Cas9"), 15)
Motif("Cas9", NNNNNNNNNNNNNNNNNNNNNGG, CCNNNNNNNNNNNNNNNNNNNNN, 21:23, 1:3, 15, true, 0)
```
"""
function setdist(motif::Motif, distance::Int)
    return Motif(motif.alias, motif.fwd, motif.rve, 
        motif.pam_loci_fwd, motif.pam_loci_rve, 
        distance, motif.extends5, motif.ambig_max)
end

# TODO add more motifs
const motif_db = Dict(
    "test" => Motif("test",
                    "NNNX",
                    "XXXG", true, true, 2, true, 0),
    "Cas9" => Motif("Cas9",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 3, true, 0),
    "Cpf1" => Motif("Cas12a",
                    "XXXXNNNNNNNNNNNNNNNNNNNN",
                    "TTTNXXXXXXXXXXXXXXXXXXXX", true, true, 3, false, 0)
    )


function Motif(alias::String; distance::Int = 3, ambig_max::Int = 0)
    motif = motif_db[alias]
    return setambig(setdist(motif, distance), ambig_max)
end


# all sequences are rotated in PAM-pattern fashion here
function appendPAM(offtarget::LongDNA{4}, motif::Motif)
    if motif.extends5
        return rve(motif.fwd[motif.pam_loci_fwd]) * offtarget
    end
    return motif.fwd[motif.pam_loci_fwd] * offtarget
end

# all sequences will be ready to search on forward strand here
# offtarget should be reversed if we are extend5 true
function appendPAM_forward(offtarget::LongDNA{4}, motif::Motif)
    if motif.extends5
        return offtarget * motif.fwd[motif.pam_loci_fwd]
    end
    return motif.fwd[motif.pam_loci_fwd] * offtarget
end
