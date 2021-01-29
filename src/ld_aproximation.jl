using BioSequences
using CSV
using DataFrames
using CRISPRofftargetHunter

## CODE

n = 20 # number of guides
glen = 20 # guide length
max_dist = 3

randoms = Set{String}()
while length(randoms) < n
   push!(randoms, String(getSeq(glen)))
end

deletion_perm = deletion_permutations(glen, max_dist, 1)
deletion_perm = [sort(vcat(x...)) for x in deletion_perm]
deletion_perm = [setdiff(collect(1:glen), x) for x in deletion_perm]
# 1    distance 0
# 20   distance 1
# 190  distance 2
# 1140 distance 3 only

# calculate real off-target table and sort them
real = zeros(Int, length(randoms), max_dist + 1)
d3 = zeros(Int, length(randoms))
for (i, r) in enumerate(randoms)
   r_del = Set([r[del] for del in deletion_perm])
   push!(d3, sum([length(x) == 17 for x in r_del]))

   for db in randoms
      if db != r
         db_del = Set([db[del] for del in deletion_perm])
         common = intersect(db_del, r_del)
         if length(common) > 0
            ld = glen - maximum([length(c) for c in common])
            real[i, ld + 1] += 1
            println("Distance ", ld, " ", r, " ", db)
         end
      end
   end
end
