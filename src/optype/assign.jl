"""
vtfinder determines and locates which torsional states each eigenstate most corresponds to
based on the sum of squares of the eigenvectors over each index within a given torsional state

svcs is the squares of the eigenvectors
jsd is the size of the rotatoinal basis
md is the size of the torsional basis
vtmax is the highest output torsional state
"""
function vtfinder(svcs,jsd::Int,md::Int,vtmax)
   ovrlp = zeros(md,size(svcs,2))
   for i in 1:md 
      ovrlp[i,:] = sum(svcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   vind = zeros(Int,size(svcs,1))
   cap = min(vtmax+4,md)
   for vi in 1:(vtmax+1)#cap 
      perm = sort(sortperm(ovrlp[vi,:], rev=true)[1:jsd])
      ovrlp[:,perm] .= 0.0
      ovrlp[vi,:] .= 0.0
      vind[perm] .= vi
   end
   return vind
end

"""
nfinder determines and locates which N states each eigenstate most corresponds to
based on the sum of squares of the eigenvectors over each index within a given N state

svcs is the squares of the eigenvectors
vind is the vt indices as from vtfinder
md is the size of the torsional basis
vtmax is the highest output torsional state
jsd is the size of the rotational basis
ns is the UnitRange of N values |J-S|:(J+S)
ni is the indices of the N levels within a J block
"""
function nfinder(svcs,vind,md,vtmax,jsd,ns,ni)
   vlist = collect(1:(vtmax+1))
   list = collect(1:size(svcs,1))[Bool.(sum(vind .∈ vlist',dims=2))[:]]
   tvcs = svcs[:,list] 
   ovrlp = zeros(length(ns),size(tvcs,2))
   for i in 1:length(ns) 
      nd = 2*ns[i]+1
      for m in 1:md
         frst = ni[i,1] + jsd*(m-1)
         last = ni[i,2] + jsd*(m-1)
         ovrlp[i,:] += transpose(sum(tvcs[frst:last, :], dims=1))
      end
   end
   nind = zeros(Int, size(svcs,1) )
   part = zeros(Int, length(list) )
   for i in 1:length(ns)
      nd = (2*ns[i]+1)*min((vtmax+1),md)
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:nd])
      nind[list[perm]] .= i
      ovrlp[:,perm] .= 0.0 
   end
   return nind
end

kperm(n::Int)::Vector{Int} = sortperm(powneg1.(collect(-n:n).+isodd(n)) .* collect(-n:n))
keperm(n::Int)::Vector{Int} = sortperm(sortperm(collect(-n:n), by=abs))[kperm(n)]
k2perm(n::Int)::Vector{Int} = sortperm(sortperm(collect(-n:n), by=abs))

"""
The goal of this function is to assign vt, N, & Kₐ by summing over vt then N blocks
   with K assigned by energetic ordering
"""
function ram36_2stg_assign(vecs,j,s,vtcalc,vtmax)
   ns,nd,ni,jsd = srprep(j,s)
   count = min(vtmax+4,vtcalc+1)
   svcs = abs.(vecs).^2
   vind = vtfinder(svcs,jsd,vtcalc+1,vtmax)
   nind = nfinder(svcs,vind,vtcalc+1,vtmax,jsd,ns,ni)
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for ng in 1:length(ns)
      nfilter = (nind .== ng)
      for v in 0:vtmax
         frst = jsd*(v) + ni[ng,1]
         last = jsd*(v) + ni[ng,2]
         part = col[nfilter .* (vind .== (v+1))]
#         part = part[keperm(ns[ng])]
         perm[frst:last] = part#col[filter][keperm(ns[ng])]
      end # v for
   end # n for
   perm = perm[perm .!= 0]
   return perm
end # function

"""
The goal of this function is to assign vt, N, & Kₐ by summing over vt then N blocks
   with K assigned by largest eigenvector elements
"""
function belgi_2stg_assign(vecs,j,s,vtcalc,vtmax)
   ns,nd,ni,jsd = srprep(j,s)
   count = min(vtmax+4,vtcalc+1)
   svcs = abs.(vecs).^2
   vind = vtfinder(svcs,jsd,vtcalc+1,vtmax)
###   if j==12.0
###      @show vind
###   end
   nind = nfinder(svcs,vind,vtcalc+1,vtmax,jsd,ns,ni)
###   if j==12.0
###      @show nind
###   end
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for ng in 1:length(ns)
      nfilter = (nind .== ng)
      for v in 0:vtmax
         frst = jsd*(v) + ni[ng,1]
         last = jsd*(v) + ni[ng,2]
         part = col[nfilter .* (vind .== (v+1))]
###         if (j==12.0)
###            @show part
###            @show allunique(part)
###            @show part[kfinder(svcs[part,part], ns[ng])]
###            @show kfinder(svcs[part,part], ns[ng])
###            @show allunique(part[kfinder(svcs[part,part], ns[ng])])
###         end
#         perm[frst:last] = part[sortperm(diag(svcs[part,part]))][keperm(ns[ng])]
         perm[frst:last] = part[kfinder(svcs[part,part], ns[ng])]
      end # v for
   end # n for
   perm = perm[perm .!= 0]
   return perm
end # function

function kfinder(block,n)
   out = zeros(Int, size(block,1))
   for i ∈ eachindex(out)
      out[i] = argmax(block[:,i])
      block[argmax(block[:,i]),:] .= 0.0
      block[:,i] .= 0.0
##      X = argmax(block[i,:])
##      out[X] = i
##      block[i,:] .= 0.0
##      block[:,X] .= 0.0
###      X = argmax(block)
###      out[X[2]] = X[1]
###      block[X[1],:] .= 0.0
###      block[:, X[2]] .= 0.0
   end
#   @show allunique(out)
   out = kbasis2store.(n,out)
#   @show allunique(out)
   return out
end
function kfinder2(block,n)
   out = zeros(Int, size(block,1))
   for i ∈ eachindex(out)
      X = argmax(block)
      out[X[2]] = X[1]
      block[X[1],:] .= 0.0
      block[:, X[2]] .= 0.0
   end
   out = kbasis2store.(n,out)
   return out
end
function kbasis2store(n::Int,x::Int)
   2*(n+1 - x)*(x < n+1) + (2*(x-n-1) + 1)*(x ≥ n+1)
end

"""
Expect-expect-overlap
-find & group ⟨ψ|ℋ_t|ψ⟩
--- this should be equal to ⟨ψ|Diag(E_ts)|ψ⟩
--- diag( sand( Diagonal( kron(E_ts,ones(ψ.R.lng)) ), vecs ) )
-find & group ⟨ψ|N²|ψ⟩ with in ℋ_t group
--- diag( sand( kron(I(vtcalc+1), N2(ψ,2,0)), vecs ) )
-sort by energy with in N²ℋ_t group

Right now the code is near exactly copied from v1 but needs to be cleaned up & checked
"""
function eeo_assign(vecs,j,s,vtcalc,vtmax)
   ns,nd,ni,jsd = srprep(j,s)
   count = min(vtmax+4,vtcalc+1)
   list = tor_expect()
   vals = vals[list]
   vecs = vecs[:,list]
   list = n_expect()
   vals = vals[list]
   vecs = vecs[:,list]
   return vals, vecs   
end

function tor_expect()
end
function n_expect()
   n2ex = diag(sand( kron(I(ctrl.vtcalc+1), N2(ψ.R,1,0)) ))
   nlist = zeros(Int,ψ.R.lng)
   for i ∈ 1:length(ψ.R.N)
      n = ψ.R.N[i]
      part = sort(sortperm(n2ex .- eh2(n), by=abs)[1:nd[i]])
      part = part[koverlap()]
      nlist[ni[i,1]:ni[i,2]] = part
      n2ex[part] .= 0.0
   end
   return nlist
end
function koverlap()
   avecs = abs.(vecs[ns:nf,part])
   piece = zeros(Int,nd)
   for j ∈ 1:nd
      a,b = Tuple(argmax(avecs))
      piece[a] = b
      avecs[:,b] *= 0.0
      avecs[a,:] *= 0.0
   end
   piece[1:(n+isodd(n))] = sort(piece[1:(n+isodd(n))],rev=true)
   piece[(n+isodd(n)+1):end] = sort(piece[(n+isodd(n)+1):end])
   return piece
end



