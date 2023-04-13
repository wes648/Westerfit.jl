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
   for v in 0:(vtmax+2) #THIS HAS TO BE SERIAL DON'T SIMD THIS ONE FUTURE WES
      mg = mcalc + vt2m(v) + 1
      perm = sort(sortperm(ovrlp[mg,:], rev=true)[1:jsd])
      mind[perm] .= mg
   end
   #println(mind)
   return mind
end
#4/9 11:30pm Something is fucked up for spin+low barrier. FUCK QNS
function mfinderv2(svcs,nind,ns,jsd,md,mcalc,vtmax)
   mind = zeros(Int,size(svcs,1))
   for ng in 1:length(ns)
      n = ns[ng]
      mind = mfinderforagivenn(svcs,mind,nind,ng,n,jsd,md,mcalc,vtmax)
   end
   return mind
end
function mfinderforagivenn(svcs,mind,nind,ng,n,jsd,md,mcalc,vtmax)
   list = collect(1:size(svcs,1))[nind .== ng]
   tvcs = svcs[:,list]
   ovrlp = zeros(md,size(tvcs,2))
   @simd for i in 1:md
      ovrlp[i,:] = sum(tvcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   part = zeros(Int,length(list))
   nd = 2*n+1
   for v in 0:vtmax
      mg = mcalc + vt2m(v) + 1
      perm = sort(sortperm(ovrlp[mg,:], rev=true)[1:nd])
      ovrlp[mg,:] .= 0.0 #prevents reassigning
      ovrlp[:,perm] .= 0.0 #prevents reassigning
      #println(perm)
      part[perm] .= mg
   end
   #println(size(part))
   #println(size(mind[list]))
   mind[list] .= part
   return mind
end

function nfinder(svcs,vtmax,md,jd,sd,ns,ni)
#this needs to be fully reworked it is only looking at the top part of the vector
   jsd = jd*sd
   ovrlp = zeros(length(ns),size(svcs,1))
   for i in 1:length(ns) 
      nd = 2*ns[i]+1
   for m in 1:md
      frst = ni[i,1] + jsd*(m-1)
      last = ni[i,2] + jsd*(m-1)
      ovrlp[i,:] += transpose(sum(svcs[frst:last, :], dims=1))
   end
   end
   nind = zeros(Int,size(svcs,1))
   #println(ovrlp)
   #ovrlp = argmax(ovrlp,dims=1)
   #println(ovrlp)
   #sortperm to grab the theoretical maximum states per n
   #sort the perm to grab the lowest vts worth of states
   #for i in 1:length(nind)
   #   nind[i] = ovrlp[i][1]
   #end
   for i in 1:length(ns)
   nd = (2*ns[i]+1)
   count = min(nd*(vtmax+3),nd*(md))
   vlimit = min(vtmax+3,md-1) 
   for v in 0:vlimit
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:count])[1:count]
      nind[perm] .= i
      #println(perm)
      ovrlp[i,:] .= 0.0
      ovrlp[:,perm] .= 0.0 
   end
   end
   #println(nind)
   return nind
end

kperm(n::Int)::Array{Int} = sortperm(Int.(cospi.(collect(-n:n).+isodd(n))) .* collect(-n:n))
keperm(n::Int)::Array{Int} = sortperm(sortperm(collect(-n:n), by=abs))[kperm(n)]

function ramassign(vecs,j::Float64,s::Float64,mcalc::Int,σt::Int,vtmax)
   svcs = abs.(vecs .* vecs)
   jd = Int(2.0*j) + 1
   sd = Int(2.0*s) + 1
   ns, nd, ni, jsd = srprep(j,s)
   #println(ns)
   #println(ni)
   md = 2*mcalc + 1 + 1*(σt==2)
   nind = nfinder(svcs,vtmax,md,jd,sd,ns,ni)
   if mcalc > 0
      mind = mfinderv2(svcs,nind,ns,jsd,md,mcalc,vtmax)
   else
      mind = ones(size(nind))
   end
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for ng in 1:length(ns)
      filter = (nind .== ng)
      for v in 0:vtmax
         mg = mcalc + vt2m(v) + 1
         filter .*= (mind .== mg)
         frst = jsd*(mg-1) + ni[ng,1]
         last = jsd*(mg-1) + ni[ng,2]
         #println("first = $frst, last = $last")
         #println(col[filter])
         #println("J = $j, ng = $ng")
         part = col[filter]
         part = part[keperm(ns[ng])]
         perm[frst:last] = part#col[filter][keperm(ns[ng])]
      end
   end
   perm = perm[perm .!= 0]
   return perm
end


### JACOBI
#Reorganize matrix to better match my conception of the quantum numbers
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

