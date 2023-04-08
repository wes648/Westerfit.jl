"""
This is where the assignment routines are located for westerfit.
   The first is a simple assigner based on RAM36. There are also the start of a 
   Jacobi eigenvalue routine based version. It is not done yet
"""

### SIMPLE
#sum across each m to find dominate torsional state
function mfinder(svcs,jsd::Int,md::Int,mcalc,vtmax)
   ovrlp = zeros(md,size(svcs,1))
   @simd for i in 1:md 
      ovrlp[i,:] = sum(svcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   #println(ovrlp)
   #ovrlp = argmax(ovrlp,dims=1)
   mind = zeros(Int,size(svcs,1))
   #@simd for i in 1:length(mind)
   #   mind[i] = ovrlp[i][1]
   #end
   for v in 0:vtmax #THIS HAS TO BE SERIAL DON'T SIMD THIS ONE FUTURE WES
      mg = mcalc + vt2m(v) + 1
      perm = sort(sortperm(ovrlp[mg,:], rev=true))
      mind[perm] .= mg
   end
   return mind
end
function nfinder(svcs,jd,sd)
   ovrlp = zeros(sd,size(svcs,1))
   @simd for i in 1:sd 
      frst = 1 + (jd-sd-1)*(i-1)
      last = (jd+1)*(i-1) + jd - sd + 1
      ovrlp[i,:] = sum(svcs[frst:last, :], dims=1)
   end
   ovrlp = argmax(ovrlp,dims=1)
   nind = zeros(Int,size(svcs,1))
   @simd for i in 1:length(nind)
      nind[i] = ovrlp[i][1]
   end
   return nind
end
keperm(n) = sortperm(sortperm(collect(-n:n), by=abs))[kperm(n)]
function ramassign(vecs,j::Float64,s::Float64,mcalc::Int,σt::Int,vtmax)
   svcs = vecs .* vecs
   jd = Int(2.0*j) + 1
   sd = Int(2.0*s) + 1
   md = 2*mcalc + 1 + 1*(σt==2)
   mind = mfinder(svcs,jd*sd,md,mcalc,vtmax)
   nind = nfinder(svcs,jd,sd)
   ns = Δlist(j,s)
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for v in 0:vtmax
      mg = mcalc + vt2m(v) + 1
      filter = (mind .== mg)
      for ng in 1:sd
         filter .*= (nind .== ng)
         frst = jd*sd*(mg-1) + 1 + (jd-sd-1)*(ng-1)
         last = jd*sd*(mg-1) + (jd+1)*(ng-1) + jd - sd + 1
         perm[frst:last] = col[filter][keperm(ns[ng])]
      end
   end
   perm = perm[perm .!= 0]
   return perm
end


### JACOBI
#Reorganize matrix to better match my conception of the quantum numbers
function kperm(n::Int)::Array{Int}
   sortperm(Int.(cospi.(collect(-n:n).+isodd(n))) .* collect(-n:n))
end
function kperm(j,s)::Array{Int}
   perm = zeros(Int,Int((2*j+1)*(2*s+1)))
   shift = 0
   for n in Δlist(j,s)
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm
end
function kperm(j,s,shift::Int,jsd::Int,Δl::Array,perm)::Array{Int}
   #perm = zeros(Int,jsd)
   for n in Δl
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm#, shift
end
function kperm(j,s,m)
   jsd = Int((2*j+1)*(2*s+1))
   Δlst = Δlist(j,s)
   shift = 0
   perm = zeros(Int,jsd*(2*m+1))
   for i in 0:(2*m)
      perm = kperm(j,s,shift,jsd,Δlst,perm)
      shift += jsd
   end
   return perm
end

#This takes the eigenvectors & builds the permutation
assignperm(vec) = sortperm([iamax(vec[:,i]) for i in 1:size(vec,2)])

