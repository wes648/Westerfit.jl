"""
This collects functions that were previously used in the code but I seem to have
replaced or otherwise abandonded. Functions that I have future plans for are
included with their related functions
"""


function qngen(n,m,σ)
   #depreciated
   nd = 2*n+1
   md = 2*m+1
   narray = fill(n,nd*md)
   karray = kron(ones(md),collect(-n:n))
   marray = kron(NFOLD .* collect(-m:m) .+ σ,ones(nd))
   σarray = fill(σ,nd*md)
   out = hcat(narray,karray,marray,σarray)
end
function qngen(j,s,m,σ)
   #depreciated
   nlist = Δlist(j,s)
   out = zeros(0,4)
   for i in 1:length(nlist)
      out = vcat(out,qngen(nlist[i],m,σ))
   end
   jlist = fill(j,size(out)[1])
   out = hcat(jlist,out)
   return out
end

function vecpadder(ns,degns,offst,nm,vtmi,vtc,vecs)
   #depreciated
   partial = Array{Float64}(undef,0,sum(degns))
   for i in 1:length(ns)
      pad = (nm - ns[i])
      zpad = zeros(Float64, pad, sum(degns))
      for v in 0:(vtc)
         temp = vecs[offst[i]+v*degns[i]+1:offst[i]+(v+1)*degns[i],:]
         temp = vcat(zpad,temp,zpad)
         partial = vcat(partial,temp)
      end
   end
   return partial
end

function greaterof(x,y)
"""
I don't remember why I needed this but it returns the greater of two input values.
"""
   if x > y
      return x
   else
      return y
   end
end

function intmat(jb,nb,jk,nk,s,k)
   sqn = length(k)
   mat = zeros(Float64,sqn,sqn)
   for x in 1:sqn
   for y in x:sqn
      @inbounds mat[x,y] = intelem(jb,nb[y],k[y],s,jk,nk[x],k[x])
   end
   end
   return Symmetric(mat)
end
function intbuild(jmax,mcalc,jb,jk,s)
   ns = Int(jmax + s)
   nbs = Δlist(jb,s)
   nks = Δlist(jk,s)
   nbarray = kron(ones((2*mcalc+1)*(2*ns+1)),nbs[1])
   nkarray = kron(ones((2*mcalc+1)*(2*ns+1)),nks[1])
   karray = kron(ones(2*mcalc+1),collect(-ns[1]:ns[1]))
   for i in 2:length(nbs)
      nbarray = vcat(nbarray,kron(ones((2*mcalc+1)*(2*ns+1)),nbs[i]))
      nkarray = vcat(nkarray,kron(ones((2*mcalc+1)*(2*ns+1)),nks[i]))
      karray = vcat(karray,kron(ones(2*mcalc+1),collect(-ns:ns)))
   end
   mat = intmat(jb,nbarray,jk,nkarray,s,karray)
   return mat
end
function intcalc(jmax,mcalc,s)
   sjmd = (2*s+1)*(2*jmax+1)*(2*mcalc+1)
   jind = convert(Int,jmax+s)
   μmat = zeros(Float64,sjmd,sjmd,jind,jind)
   for x in 1:jind
   for y in x:jind
      μmat[:,:,x,y] = intbuild(jmax,y-s,x-s,s)
      μmat[:,:,y,x] = μmat[:,:,x,y]
   end
   end
   return μmat
end
