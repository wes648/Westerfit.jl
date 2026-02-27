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
jd is the degeneracy of ψ.R.J
sd is the degeneracy of ψ.R.S
ns is the UnitRange of N values |J-S|:(J+S)
ni is the indices of the N levels within a J block
"""
function nfinder(svcs,vind,md,vtmax,jd,sd,ns,ni)
   jsd = jd*sd
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
   nind = zeros(Int,size(svcs,1))
   part = zeros(Int,length(list))
   for i in 1:length(ns)
      nd = (2*ns[i]+1)*min((vtmax+1),md
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:nd])
      nind[list[perm]] .= i
      ovrlp[:,perm] .= 0.0 
   end
   return nind
end

"""
The goal of this function is to assign vt, N, & Kₐ by summing over vt and N blocks
"""
function ram36_2stg_assign(vecs,j,s,vtcalc,vtmax)
   ns,nd,ni,jsd = srprep(j,s)
   count = min(vtmax+4,vtcalc+1)
   svcs = abs.(vecs).^2
   vind = vtfinder(svc,jsd,)
end

