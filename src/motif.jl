"""
`Motif(
    alias::String,
    fwdmotif::String, 
    fwdpam::String, 
    forward_strand::Bool = true, 
    reverse_strand::Bool = true, 
    distance::Int = 4, 
    extends5::Bool = true,
    ambig_max::Int = 5)`

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

Alignments will be performed from opposite to the extension direction (which is deifned by extend5).

# Examples
```julia-repl
Motif('Cas9')
Motif('Cas9', 'NNNNNNNNNNNNNNNNNNNNXXX', 'XXXXXXXXXXXXXXXXXXXXNGG'. true, true, 3, true, 5)
```
"""
struct Motif
    alias::String
    fwd::LongDNASeq
    rve::LongDNASeq
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
function removepam(seq::LongDNASeq, pam::UnitRange{<:Integer})
    x = copy(seq)
    deleteat!(x, pam)
    return x
end


function Motif(alias::String,
    fwdmotif::String, fwdpam::String,
    forward_strand::Bool = true, reverse_strand::Bool = true,
    distance::Int = 4, extends5::Bool = true, ambig_max::Int = 4)
    if length(fwdmotif) != length(fwdpam)
        throw("fwd_motif and fwd_pam have to have the same length!")
    end
    merge = combinestrings(fwdmotif, fwdpam)

    if forward_strand
        # where is PAM located?
        pam_loci_fwd = findall(r"[^X]+", fwdpam)[1]
        fwd = LongDNASeq(merge)
    else
        pam_loci_fwd = UnitRange{Int64}()
        fwd = LongDNASeq("")
    end

    if reverse_strand
        pam_loci_rve = findall(r"[^X]+", reverse(fwdpam))[1]
        rve = reverse_complement(LongDNASeq(merge))
    else
        pam_loci_rve = UnitRange{Int64}()
        rve = LongDNASeq("")
    end

    return Motif(alias, fwd, rve, pam_loci_fwd, pam_loci_rve, distance, extends5, ambig_max)
end


"
Calculate what is the length of the motif, with extension, but without PAM.
Effectively, size of the off-target on the genome.
"
function length_noPAM(motif::Motif)
    fwd_len = length(motif.fwd) - length(motif.pam_loci_fwd)
    rve_len = length(motif.rve) - length(motif.pam_loci_rve)
    return max(fwd_len, rve_len)
end


function length(motif::Motif)
    return max(length(motif.fwd), length(motif.rve))
end


function setambig(motif::Motif, ambig::Int)
    return Motif(motif.alias, motif.fwd, motif.rve, 
        motif.pam_loci_fwd, motif.pam_loci_rve, 
        motif.distance, motif.extends5, ambig)
end


function setdist(motif::Motif, distance::Int)
    return Motif(motif.alias, motif.fwd, motif.rve, 
        motif.pam_loci_fwd, motif.pam_loci_rve, 
        distance, motif.extends5, motif.ambig_max)
end



# function print_rgb(io::IO, t::String, r = 235, g = 79, b = 52)
#     print(io, "\e[1m\e[38;2;$r;$g;$b;249m", t)
# end
#
# function print_pam(io::IO, s::LongSequence, pam_idx::Vector{Int})
#     for i in eachindex(s)
#         if i in pam_idx
#             print_rgb(io, s[i])
#         else
#             print_rgb(io, s[i], 52, 207, 235)
#         end
#     end
# end
#
# function Base.show(io::IO, mime::MIME"text/plain", motif::Motif)
#     compact = get(io, :compact, false)
#
#     if !compact
#         alias = join(["\nMotif:", motif.alias, "\n"])
#         print_rgb(io, alias, 176, 176, 176)
#         print_rgb(io, "\nForward pattern:\n", 176, 176, 176)
#         pam_idx = vcat(map(collect, motif.pam_loci_fwd)...)
#         print_rgb(io, string(motif.fwd.seq), pam_idx)
#         print_rgb(io, "\nReverse pattern:\n", 176, 176, 176)
#         pam_idx = vcat(map(collect, motif.pam_loci_rve)...)
#         print_rgb(io, string(motif.rve.seq), pam_idx)
#     else
#         show(io, motif)
#     end
# end

# TODO add more motifs
const motif_db = Dict(
    "Cas9" => Motif("Cas9",
                    "NNNNNNNNNNNNNNNNNNNNXXX",
                    "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 4, true, 0),
    "Cpf1" => Motif("Cas12a",
                    "XXXXNNNNNNNNNNNNNNNNNNNN",
                    "TTTNXXXXXXXXXXXXXXXXXXXX", true, true, 4, false, 0)
    )


function Motif(alias::String; distance::Int = 4, ambig_max::Int = 0)
    motif = motif_db[alias]
    return setambig(setdist(motif, distance), ambig_max)
end
