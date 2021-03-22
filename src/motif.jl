"
Motif defines what we search on the genome,
what can be identified as an off-target.

`distance` defines how many bases extra we need on the
off-target sequence for the alignment.

For example for Cas9- 20bp-NGG with up to 4 mm
we need to have 24bp-NGG patterns on the gneome.
"
struct Motif
    alias::String
    fwd::LongDNASeq
    rve::LongDNASeq
    pam_loci_fwd::Vector{UnitRange{<:Integer}}
    pam_loci_rve::Vector{UnitRange{<:Integer}}
    distance::Int
    extends5::Bool
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

function combinestrings(rule::Function, s1::String, s2::String)
    if (length(s1) != length(s2))
        error("Unequal lengths.")
    end
    return join([rule(s1[i], s2[i]) for i in eachindex(s1)])
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


"
Constructor for Motif that will be found in the reference.

`alias`    - alias of the motif for easier identification e.g. Cas9
`fwdmotif` - Motif that indicates where is PAM inside `fwdpam`.
             For example for Cas9 it is 20*N + XXX:
             NNNNNNNNNNNNNNNNNNNNXXX
`fwdpam`   - Motif in 5'-3' that will be matched on the reference (without the X).
             For example for Cas9 it is 20*N + NGG:
             XXXXXXXXXXXXXXXXXXXXNGG
`forward`  - If false will not match to the forward reference strand.
`reverse`  - If false will not match to the reverse reference strand.
`distance` - How many extra nucleotides are needed for a search? This
             will indicate within what distance we can search for off-targets.
`extend5`  - Defines how off-targets will be aligned to the guides and where
             extra nucleotides will be added for alignment within distance. Whether
             to extend in the 5' and 3' direction.

Example for Cas9 where we want to search for off-targets within distance of 4:
alias:    Cas9
fwdmotif:     NNNNNNNNNNNNNNNNNNNNXXX
fwdpam:       XXXXXXXXXXXXXXXXXXXXNGG
forward:  true
reverse:  true
distance: 4
extend5:  true

This will search these motifs on the genome (forward strand)
           NNNNNNNNNNNNNNNNNNNNNNNNNGG
           EEEE

           CCNNNNNNNNNNNNNNNNNNNNNNNNN      (reverse strand)
                                  EEEE
and treat those 4 E nucleotides as an extension, only used for alignment
purposes, and alignments will start from opposite to the `extend5` direction.
"
function Motif(alias::String,
    fwdmotif::String, fwdpam::String,
    forward_strand::Bool = true, reverse_strand::Bool = true,
    distance::Int = 4, extends5::Bool = true)
    if length(fwdmotif) != length(fwdpam)
        throw("fwd_motif and fwd_pam have to have the same length!")
    end
    if extends5
        fwdmotif = repeat("N", distance) * fwdmotif
        fwdpam = repeat("X", distance) * fwdpam
    else
        fwdmotif = fwdmotif * repeat("N", distance)
        fwdpam = fwdpam * repeat("X", distance)
    end
    merge = combinestrings(notX, fwdmotif, fwdpam)

    if forward_strand
        # where is PAM located?
        pam_loci_fwd = findall(r"[^X]+", fwdpam)
        fwd = LongDNASeq(merge)
    else
        pam_loci_fwd = Vector{UnitRange{Int64}}()
        fwd = LongDNASeq("")
    end

    if reverse_strand
        pam_loci_rve = findall(r"[^X]+", reverse(fwdpam))
        rve = reverse_complement(LongDNASeq(merge))
    else
        pam_loci_rve = Vector{UnitRange{Int64}}()
        rve = LongDNASeq("")
    end

    return Motif(alias, fwd, rve, pam_loci_fwd, pam_loci_rve, distance, extends5)
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
                    "XXXXXXXXXXXXXXXXXXXXNGG", true, true, 4, true),
    "Cpf1" => Motif("Cas12a",
                    "XXXXNNNNNNNNNNNNNNNNNNNNNNNN",
                    "TTTNXXXXXXXXXXXXXXXXXXXXXXXX", true, true, 4, false)
    )

function Motif(alias::String)
    return motif_db[alias]
end
