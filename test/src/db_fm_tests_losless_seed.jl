#=
using CRISPRofftargetHunter
using BioSequences
using CSV
using DataFrames
using FASTX
using StatsBase
using StaticArrays

cd("test")

genome = joinpath(dirname(pathof(CRISPRofftargetHunter)), "..", 
    "test", "sample_data", "genome", "semirandom.fa")
#genome = "/home/ai/Projects/uib/crispr/chopchop_genomes/hg38v34.fa"
guides_s = Set(readlines("./sample_data/crispritz_results/guides.txt"))
guides = LongDNASeq.(guides_s)

motif = Motif("Cas9")
motif = CRISPRofftargetHunter.setdist(motif, 3)

fmidir = tempname()
mkpath(fmidir)

ldb_path = joinpath(fmidir, "linearDB")
mkpath(ldb_path)
build_linearDB("samirandom", genome, motif, ldb_path, 7)
ldb_res = search_linearDB(ldb_path, guides, 3)

fmidbpath = build_fmiDB(genome, fmidir)
pamdbpath = joinpath(fmidbpath, "pamDB.bin")
build_pamDB(fmidir, motif; storagedir = pamdbpath)

fmidbdir = fmidbpath 
genomepath = genome
detail = ""
distance = 3

pamDB = CRISPRofftargetHunter.load(pamdbpath)
if distance > pamDB.motif.distance
    throw("DIstance is too large for selected MotifPathTemplates. Max distance is" * 
        string(pamDB.motif.distance))
end

# 01*0 seed always exists when pattern is partitioned in at leasst distance + 2 parts
adj_seed_size = Int(floor(length(guides[1])/(distance + 2)))
parts_number = Int(floor(length(guides[1])/adj_seed_size))

gi = CRISPRofftargetHunter.load(joinpath(fmidbdir, "genomeInfo.bin"))
is_fa = CRISPRofftargetHunter.is_fasta(genomepath)

guides_ = copy(guides)
# reverse guides so that PAM is always on the left
# TODO CPf1 style...
if pamDB.motif.extends5
    guides_ = reverse.(guides_)
end
guides_skipmers = Base.map(x -> unique(CRISPRofftargetHunter.as_skipkmers(x, adj_seed_size)), guides_)

res = zeros(Int, length(guides_), distance + 1)
ref = open(genomepath, "r")
reader = is_fa ? FASTA.Reader(ref, index = genomepath * ".fai") : TwoBit.Reader(ref)


offt_len = CRISPRofftargetHunter.length_noPAM(pamDB.motif) + distance
# 1       23
# hit---hitPAM
# |     |  
scan_dist = offt_len - adj_seed_size
# worst case
# hithit---PAM
#       |  | = 23 - adj_seed_size * 2 = 15
# hit-hit--PAM
# |   |
# diff
#        | | 
# 23 - diff - adj_seed_size 
# adj_seed_size is the difference here 
# hit---hitPAM
#          | = 0
for (ic, chrom) in enumerate(gi.chrom)
    fmi = CRISPRofftargetHunter.load(joinpath(fmidbdir, chrom * ".bin"))
    seq = CRISPRofftargetHunter.getchromseq(is_fa, reader[chrom])
    pam_loc_fwd = pamDB.pam_loc_fwd[ic]
    pam_loc_fwd_len = length(pam_loc_fwd)
    pam_loc_rve = pamDB.pam_loc_rve[ic]
        
    for (i, gs) in enumerate(guides_skipmers)
        # FORWARD
        gs_fwd = reverse.(gs)
        hits = mapreduce(x -> CRISPRofftargetHunter.locateall(x, fmi), vcat, gs_fwd)
        hits = unique(hits)
        sort!(hits)
        diff_hits = diff(hits)
        two_are_close = diff_hits .<= scan_dist
        diff_hits = diff_hits[two_are_close]
        pushfirst!(two_are_close, 0) # we take the most right element
        hits = hits[two_are_close]
        # hitsPAM
        # |-> |
        hits .= hits .+ adj_seed_size
        check_offtarget_pam = zeros(Int, length(hits))
        check_offtarget = falses(length(hits))
        for (ih, h) in enumerate(hits)
            idx = searchsortedfirst(pam_loc_fwd, h) # find first >= h
            if idx <= pam_loc_fwd_len
                # |    |-PAM
                # 14 462 469
                if (pam_loc_fwd[idx] - h) <= (scan_dist - diff_hits[ih])
                    check_offtarget[ih] = true
                    check_offtarget_pam[ih] = pam_loc_fwd[idx]
                end
            end
        end
        check_offtarget_pam = unique(check_offtarget_pam[check_offtarget])
        check_offtarget_pam .= check_offtarget_pam .- 1
        offt_seq = map(x -> seq[(x - offt_len + 1):x], check_offtarget_pam)
        if pamDB.motif.extends5
            offt_seq = reverse.(offt_seq)
        end
        align_score = map(x -> CRISPRofftargetHunter.levenshtein(guides_[i], x, distance), offt_seq)
        for d in 0:distance
            res[i, d + 1] += sum(align_score .== d)
        end

        # REVERSE
        gs_rve = complement.(gs)
        hits = mapreduce(x -> CRISPRofftargetHunter.locateall(x, fmi), vcat, gs_rve)
        hits = unique(hits)
        sort!(hits)
        diff_hits = diff(hits)
        two_are_close = diff_hits .<= scan_dist
        diff_hits = diff_hits[two_are_close]
        push!(two_are_close, 0) # we take the most left element
        hits = hits[two_are_close]
        hits .= hits .- 1

        check_offtarget_pam = zeros(Int, length(hits))
        check_offtarget = falses(length(hits))
        for (ih, h) in enumerate(hits)
            idx = searchsortedlast(pam_loc_rve, h) # find last value <= h
            if idx != 0 
                if (h - pam_loc_rve[idx]) <= (scan_dist - diff_hits[ih])
                    check_offtarget[ih] = true
                    check_offtarget_pam[ih] = pam_loc_rve[idx]
                end
            end
        end
        check_offtarget_pam = unique(check_offtarget_pam[check_offtarget])
        check_offtarget_pam .= check_offtarget_pam .+ 1
        offt_seq = map(x -> seq[x:(x + offt_len - 1)], check_offtarget_pam)
        align_score = map(x -> CRISPRofftargetHunter.levenshtein(guides_[i], complement(x), distance), offt_seq)
        for d in 0:distance
            res[i, d + 1] += sum(align_score .== d)
        end
    end
end
close(ref)

res = DataFrame(res, :auto)
col_d = [Symbol("D$i") for i in 0:distance]
rename!(res, col_d)
res.guide = guides
sort!(res, col_d)
=#