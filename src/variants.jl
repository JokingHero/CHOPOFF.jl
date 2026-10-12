# Sequence-resolved variants for prefixHashScan. These helpers are internal;
# the public interface is search_prefixHashScan(...; vcf=...).
include("variants/vcf_read.jl")
include("variants/fragments.jl")
include("variants/enumerate.jl")
include("variants/search.jl")
