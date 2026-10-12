using Random

vtest_alternative(base) = base == 'A' ? "C" : "A"
vtest_rc(sequence) = string(reverse_complement(LongDNA{4}(sequence)))

function vtest_spacing_case(name, distance, spacing, reverse_strand)
    guide = name == "Cas9" ? VCF_TEST_GUIDE : VCF_TEST_GUIDE * "A"
    site = name == "Cas9" ? guide * "AGG" : "TTTA" * guide
    reference = "C"^40 * site * "C"^80
    a = name == "Cas9" ? 41 : 45
    b = a + spacing
    alt_a = vtest_alternative(reference[a])
    # Preserve constrained PAM bases when the second SNP reaches the PAM.
    alt_b = name == "Cas9" && b in (62,63) ? reference[b:b] : vtest_alternative(reference[b])
    edits = [vtest_edit(reference,a,a+1,alt_a,"a"),vtest_edit(reference,b,b+1,alt_b,"b")]
    first_guide = alt_a * guide[2:end]
    guides = [LongDNA{4}(first_guide)]
    if spacing < length(guide)
        both = collect(first_guide)
        both[spacing+1] = only(alt_b)
        push!(guides,LongDNA{4}(String(both)))
    end
    if reverse_strand
        n = length(reference)
        reference = vtest_rc(reference)
        edits = [vtest_edit(reference,n-v.stop+2,n-v.first+2,vtest_rc(v.alt),v.id) for v in edits]
    end
    return reference,edits,unique(guides),Motif(name;distance)
end

function vtest_random_case(rng, i)
    name = isodd(i) ? "Cas9" : "Cas12a"
    reverse_strand = isodd(div(i-1,2))
    distance = (i-1) % 5
    g = String(rand(rng,collect("ACGT"),name == "Cas9" ? 20 : 21))
    site = name == "Cas9" ? g*"AGG" : "TTTA"*g
    reference = String(rand(rng,collect("ACGT"),40)) * site * String(rand(rng,collect("ACGT"),40))
    p = name == "Cas9" ? 45 : 49
    if reverse_strand
        p = length(reference)-p+1
        reference = vtest_rc(reference)
    end
    alt = vtest_alternative(reference[p])
    guide = collect(g)
    guide[5] = only(reverse_strand ? vtest_rc(alt) : alt)
    edits = [vtest_edit(reference,p,p+1,alt,"anchor")]
    n = rand(rng,1:8)
    for j in 2:n
        q = rand(rng,35:(length(reference)-35))
        kind = rand(rng,1:5)
        len = rand(rng,1:3)
        if kind == 1
            push!(edits,vtest_edit(reference,q,q+1,vtest_alternative(reference[q]),"v$j"))
        elseif kind == 2
            insertion = String(rand(rng,collect("ACGT"),len-1)) * vtest_alternative(reference[q-1])
            push!(edits,vtest_edit(reference,q,q,insertion,"v$j"))
        elseif kind == 3
            # Canonical deletion padding differs from the last deleted base.
            while reference[q-1] == reference[q+len-1]
                q = rand(rng,35:(length(reference)-35))
            end
            push!(edits,vtest_edit(reference,q,q+len,"","v$j"))
        else
            ref = reference[q:q+len-1]
            replacement = join(vtest_alternative(base) for base in ref)
            kind == 5 && (replacement = len == 1 ? replacement^2 : replacement[1:1] * "C" * replacement[2:end])
            push!(edits,vtest_edit(reference,q,q+len,replacement,"v$j"))
        end
    end
    return reference,edits,[LongDNA{4}(String(guide))],Motif(name;distance),distance
end
