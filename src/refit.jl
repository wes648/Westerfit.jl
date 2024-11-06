
"""
This file is for re-assigning QNs from other programs by finding the
   permuation of states that gives the lowest initial rms.

Need a better idea that doesn't involve insade numbers of factorials
"""

using Combinatorics

function nlinds(n)
   snd = convert(Int, sum(2 .* collect(0:(n-1)) .+ 1)) +1
   fnd = convert(Int, sum(2 .* collect(0:n) .+ 1))
   return collect(snd:fnd)
end

function limsubtract(a,b)
   r,c, = findnz(a)
   f = sparse(r,c,ones(length(c)))
   out = a - (b .* f)
   return out
end

permlength(nmax)::Int = (nmax+1)^2

function fullpermcount(nmax)
   list = factorial.(2 .* collect(0:nmax) .+ 1)
   list = prod(list)
   return list
end

function build_nthperm(nmax,k,l)
   perm = ones(Int,l)
   @simd for n in 1:nmax
      nlist = nlinds(n)
   @inbounds perm[nlist] = nthperm(nlist,mod1(k,factorial(length(nlist))))
   end
   return perm
end
function build_nthperm(nmax,k)
   perm = build_nthperm(nmax,k,permlength(nmax))
   return perm
end

function build_allperms(nmax)
   list = fullpermcount(nmax)
   leng = permlength(nmax)
   allperms = zeros(Int,leng,list)
   @threads for k in 1:list 
      allperms[:,k] = build_nthperm(nmax,k,leng)
   end
   return allperms
end

