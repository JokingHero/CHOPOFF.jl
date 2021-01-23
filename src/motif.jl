"
we need to search using pattern with motif
but add to the structures only bases without PAMs

how to define PAM nicely?
Input is simple:

fwd style pattern with XXX where PAM
NNNNXXX
XXXXNGG - PAM
reverse strand - t/f

what we need is:
transformer seqeunce fwd -> no pam seq fwd
transformer          rve -> no pam rve

two methods bit level or DNA level removal
"
struct Motif
    alias::String
    fwd::ExactSearchQuery{LongSequence{DNAAlphabet{4}}}
    rve::ExactSearchQuery{LongSequence{DNAAlphabet{4}}}
    pam_loci_fwd::Vector{UnitRange{<:Integer}}
    pam_loci_rve::Vector{UnitRange{<:Integer}}
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

import Base.mergewith
function mergewith(rule::Function, s1::String, s2::String)
    if (length(s1) != length(s2))
        throw("Unequal lengths.")
    end
    return join([rule(s1[i], s2[i]) for i in eachindex(s1)])
end

"
Constructor for Motif that will be found in the reference.

`alias` - alias of the motif for easier identification e.g. Cas9
`fwdmotif` - Motif that indicates where is PAM inside `fwdpam`.
             For example for Cas9 it is 20*N + XXX:
             NNNNNNNNNNNNNNNNNNNNXXX
`fwdpam` - Motif in 5'-3' that will be matched on the reference (without the X).
           For example for Cas9 it is 20*N + NGG:
             XXXXXXXXXXXXXXXXXXXXNGG
`forward` - If false will not match to the forward reference strand.
`reverse` - If false will not match to the reverse reference strand.
"
function Motif(alias::String,
    fwdmotif::String, fwdpam::String,
    forward_strand = true, reverse_strand = true)
    if length(fwdmotif) != length(fwdpam)
        throw("fwd_motif and fwd_pam have to have the same length!")
    end
    merge = mergewith(notX, fwdmotif, fwdpam)
    merge = LongDNASeq(merge)

    if forward_strand
        # where is PAM located?
        pam_loci_fwd = findall(r"[^X]+", fwdpam)
        fwd = ExactSearchQuery(merge)
    else
        pam_loci_fwd = Vector{UnitRange{Int64}}()
        fwd = ExactSearchQuery(LongDNASeq(""))
    end

    if reverse_strand
        pam_loci_rve = findall(r"[^X]+", reverse(fwdpam))
        rve = ExactSearchQuery(reverse_complement(merge))
    else
        pam_loci_rve = Vector{UnitRange{Int64}}()
        rve = ExactSearchQuery(LongDNASeq(""))
    end

    return Motif(alias, fwd, rve, pam_loci_fwd, pam_loci_rve)
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

# motif database - check this with other software
const motif_db = Dict(
    "Cas9" => Motif("Cas9", "NNNNNNNNNNNNNNNNNNNNXXX", "XXXXXXXXXXXXXXXXXXXXNGG"),
    "Cpf1" => Motif("Cpf1", "XXXXNNNNNNNNNNNNNNNNNNNNNNNN", "TTTNXXXXXXXXXXXXXXXXXXXXXXXX"),
    "CasX" => Motif("CasX", "XXXXNNNNNNNNNNNNNNNNNNNNNNNN", "TTCNXXXXXXXXXXXXXXXXXXXXXXXX"),
    "Cas13" => Motif("Cas13", "XNNNNNNNNNNNNNNNNNNNNNNNNNNN", "HXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    )

function Motif(alias::String)
    return motif_db[alias]
end
