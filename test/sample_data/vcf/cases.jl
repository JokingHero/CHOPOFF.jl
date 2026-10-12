# Explicit half-open edits feed the fixture writer and independent haplotype
# oracle. No production VCF parser or coordinate helper participates here.
const VCF_TEST_GUIDE = "ACGTGTCAGTACGATCGTAC"
const VCF_TEST_REFERENCE = "T"^40 * VCF_TEST_GUIDE * "AGG" * "T"^80
const VCF_TEST_IUPAC = Dict('R'=>"AG", 'Y'=>"CT", 'S'=>"CG", 'W'=>"AT", 'K'=>"GT",
    'M'=>"AC", 'B'=>"CGT", 'D'=>"AGT", 'H'=>"ACT", 'V'=>"ACG", 'N'=>"ACGT")

function vfixture(id, purpose, reference, edits; guide = VCF_TEST_GUIDE, motif = "Cas9",
    distance = 2, witness = [first(edits).id], variant = true, chrom = "chr1", requirement = "V04")
    return (; id, purpose, reference, edits, guides = [LongDNA{4}(guide)], motif,
        distance, witness, variant, chrom, requirement)
end

function vfixture_motif(case; distance = case.distance, ambig_max = 0)
    case.motif == "pamless" && return Motif("pamless", "N"^20, "X"^20, true, true, distance, true, ambig_max)
    case.motif == "internal" && return Motif("internal", "N"^10 * "XXX" * "N"^10,
        "X"^10 * "AGG" * "X"^10, true, true, distance, false, ambig_max)
    case.motif == "forward_only" && return Motif("forward_only", "N"^20 * "XXX",
        "X"^20 * "AGG", true, false, distance, true, ambig_max)
    return Motif(case.motif; distance, ambig_max)
end

function vfixture_cases()
    r, g = VCF_TEST_REFERENCE, VCF_TEST_GUIDE
    e(a, b, alt, id; kwargs...) = vtest_edit(r, a, b, alt, id; kwargs...)
    cases = [
        vfixture("mixed_choices", "REF/ALT combinations of SNP, insertion and deletion", r,
            [e(45,46,"A","snp"), e(51,51,"C","ins"), e(54,55,"","del")]; witness=["snp","ins"]),
        vfixture("mnp_replacement", "MNP and unequal-length replacement in one site", r,
            [e(46,48,"GA","mnp"), e(53,54,"CA","replace")]; distance=4, witness=["mnp","replace"]),
        vfixture("nontransitive_conflicts", "An overlapping deletion does not merge independent choices", r,
            [e(44,48,"","del1"), e(45,46,"A","inside"), e(47,50,"","del2")]; distance=4, witness=["inside","del2"], requirement="V05"),
        vfixture("same_gap", "Two independent insertions at one gap are exclusive", r,
            [e(46,46,"C","gapC"), e(46,46,"A","gapA")]; requirement="V05"),
        vfixture("gap_boundaries", "Insertions at both boundaries of a deletion remain compatible", r,
            [e(46,48,"","del"), e(46,46,"A","left"), e(48,48,"T","right")]; distance=4, witness=["del","left","right"], requirement="V05"),
        vfixture("adjacent_deletions", "Adjacent deletion junctions remain distinct alleles", r,
            [e(47,48,"","left"), e(48,49,"","right")]; witness=["left","right"], requirement="V05"),
        vfixture("siblings", "SNP and insertion ALTs from one record cannot combine", r,
            [e(45,46,"A","multi";group="source"), e(46,46,"C","multi";group="source")]; requirement="V05"),
        vfixture("split_sources", "Independent split records permit the SNP-plus-insertion walk", r,
            [e(45,46,"A","multi";group="snp"), e(46,46,"C","multi";group="ins")]; witness=["multi","multi"], requirement="V05"),
        vfixture("restore_reference", "Distinct edits reconstruct reference and retain both annotations", r,
            [e(47,48,"","remove"), e(47,47,"C","restore")]; distance=0, witness=["remove","restore"]),
        vfixture("noop", "An explicit reference-equivalent ALT remains annotated", r,
            [e(45,46,"G","same")]; distance=0),
        vfixture("missing_id", "Missing ID falls back to allele coordinates", r,
            [e(45,46,"A",".")]; witness=["chr1:45:G>A"], requirement="V10"),
        vfixture("pam_destroy", "PAM destruction leaves only the annotated reference hit", r,
            [e(62,63,"A","pam")]; distance=0, variant=false, witness=String[], requirement="V08"),
    ]
    created = r[1:61] * "A" * r[63:end]
    push!(cases, vfixture("pam_create", "ALT creates a PAM absent from the reference", created,
        [vtest_edit(created,62,63,"G","pam")]; distance=0, requirement="V08"))
    for (name, replacement) in (("pam_guide_mnp", "AT"), ("pam_guide_replacement", "ACT"))
        push!(cases, vfixture(name, "Replacement crosses the guide/PAM boundary", r,
            [e(60,62,replacement,"junction")]; distance=3, requirement="V08"))
    end
    for (name, inserted) in (("long_insertion", g*"AGG"*"T"^35*g*"AGG"),
            ("long_insertion_reverse", "C"^70*string(reverse_complement(LongDNA{4}(g*"AGG")))*"C"^75),
            ("long_insertion_both", (g*"AGG"*"C"^31)^4*string(reverse_complement(LongDNA{4}(g*"AGG")))*"C"))
        push!(cases, vfixture(name, "Sites inside an insertion retain distinct origin offsets", r,
            [e(90,90,inserted,"long")]; distance=0, requirement="V09"))
    end
    for (name, first, stop, alt) in (("endpoint_insertion",61,61,"AGG"),
            ("guide_deletion",49,51,""), ("guide_insertion",49,49,"T"))
        push!(cases, vfixture(name, "Indel at an alignment boundary exercises bulges and endpoints", r,
            [e(first,stop,alt,"boundary")]; distance=3, requirement="V09"))
    end
    distant = "T"^40 * g[1:10] * "C"^70 * g[11:20] * "AGG" * "T"^40
    push!(cases, vfixture("deletion_bridge", "A deletion brings a distant SNP into the same site", distant,
        [vtest_edit(distant,46,47,"A","snp"),vtest_edit(distant,51,121,"","bridge")]; witness=["snp","bridge"], requirement="V07"))
    edge = g * "AGG" * "T"^50
    push!(cases, vfixture("left_edge", "First-base insertion/deletion and nearby SNP", edge,
        [vtest_edit(edge,1,1,"C","before"),vtest_edit(edge,1,2,"","delete"),vtest_edit(edge,3,4,"A","near")]; witness=["near"], requirement="V07"))
    edge = "T"^50 * g * "AGG"
    push!(cases, vfixture("right_edge", "An allele touches the last chromosome base", edge,
        [vtest_edit(edge,length(edge),length(edge)+1,"G","last")]; distance=0, requirement="V07"))
    rev = string(reverse_complement(LongDNA{4}(r)))
    push!(cases, vfixture("reverse_snp", "Reverse-strand allele and endpoint mapping", rev,
        [vtest_edit(rev,99,100,"G","reverse")]; requirement="V08"))
    for (name, motif, guide, site, p) in (("cas12a","Cas12a",g*"A","TTTA"*g*"A",49),
            ("generic_nga","Cas9_NGA",g,g*"AGA",45),
            ("pamless","pamless",g,g,45),
            ("internal_pam","internal",g,g[1:10]*"AGG"*g[11:20],45),
            ("forward_only","forward_only",g,g*"AGG",45))
        ref = "C"^40 * site * "C"^60
        alt = ref[p] == 'A' ? "C" : "A"
        push!(cases, vfixture(name, "Motif-specific coordinates and ALT scanning", ref,
            [vtest_edit(ref,p,p+1,alt,"motif")]; guide, motif, requirement="V08"))
    end
    for code in sort(collect(keys(VCF_TEST_IUPAC)))
        edits = [e(45,46,string(base),"amb";original=string(code),group="amb") for base in VCF_TEST_IUPAC[code]]
        push!(cases,vfixture("iupac_$code", "Every concrete ALT realization retains its identity",r,edits;requirement="V06"))
    end
    edits = [e(45,47,string(a,b),"amb2";original="RY",group="amb2") for a in "AG" for b in "CT"]
    push!(cases,vfixture("iupac_two_positions","Two-position IUPAC expansion stays allele-specific",r,edits;requirement="V06"))
    return cases
end

function vfixture_invalid()
    row(pos, ref, alt; chrom="chr1") = (chrom,pos,"bad",ref,alt,".","PASS",".")
    return [
        (id="unknown_contig", records=[row(45,"G","A";chrom="absent")], message="absent from the reference"),
        (id="ref_mismatch", records=[row(45,"T","A")], message="REF mismatch"),
        (id="outside_reference", records=[row(144,"A","C")], message="outside the reference"),
        (id="zero_position", records=[row(0,"A","C")], message="outside the reference"),
        (id="ref_crosses_end", records=[row(143,"TT","C")], message="outside the reference"),
        (id="unsorted", records=[row(46,"T","C"),row(45,"G","A")], message="not sorted"),
        (id="repeated_contig", records=[row(45,"G","A"),row(45,"G","A";chrom="chr2"),row(46,"T","C")], message="not contiguous"),
        [(id="unsupported_$name", records=[row(45,"G",alt)], message="Unsupported VCF ALT")
            for (name,alt) in (("del","<DEL>"),("dup","<DUP>"),("inv","<INV>"),
                ("ins","<INS>"),("cnv","<CNV>"),("bnd","G]chr2:45]"),("letter","Z"))]...,
    ]
end

function vfixture_contracts()
    annotations = (
        id="annotations",
        header=["##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT frequency\">",
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
            "##INFO=<ID=PAIR,Number=2,Type=Integer,Description=\"Pair\">",
            "##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"Flag\">",
            "##INFO=<ID=NOTE,Number=.,Type=String,Description=\"Text\">",
            "##INFO=<ID=MISSING,Number=1,Type=String,Description=\"Absent\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"ALT fraction\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">"],
        samples=["sample1","sample%2"],
        records=[("chr1",45,"id%q;id2","G","A,C",".","PASS",
            "AF=0.1,0.2;AD=100,10,20;PAIR=7,8;FLAG;NOTE=a,b\"c%|","GT:AF:AD:PL",
            "0|1:0.25,0.75:8,2,6:0,1,2,3,4,5","./.:.:.:.")])
    markers = (id="ignored_markers",header=String[],samples=String[],records=[
        ("chr1",45,"no_alt","G",".",".","PASS","."),
        ("chr1",46,"spanning","T","*",".","PASS","."),
        ("chr1",47,"mixed","C","<NON_REF>,A",".","PASS","."),
        ("chr1",48,"nonref","A","<*>",".","PASS",".")])
    return [annotations, markers]
end
