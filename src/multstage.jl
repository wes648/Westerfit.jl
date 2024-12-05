"""
This file is for testing variations on the theme of multiple torsional stages
"""

#include("/home/wes/files/westerfit/src/main.jl")

function build_h1(n,mc,σ)
   f = 5
   v = 200
   ms = collect(-(3mc+σ):3:(3mc+σ))
   ht = spdiagm(0=>(f.* ms.^2 .+ 0.5*v), 1=>fill(-0.25v,length(ms)-1), -1=>(fill(-0.25v,length(ms)-1))) .* csl
   hr = hrotest([1.75;1.25;0.125;.08],n)
   rho = -2*f*csl*0.5 .* kron(Diagonal(ms),Diagonal(collect(-n:n)))
   out = dropzeros!(kron(ht,I(2n+1)) + kron(I(length(ms)),hr) + rho)
   out = kron(ur(mc),ur(n)) * out * kron(ur(mc),ur(n))
   return dropzeros!(out) ./csl
end

function build_h2s(n,mc,mm,σ)
   f = 5
   v = 200
   ms = collect(-(3mc+σ):3:(3mc+σ))
   ht = spdiagm(0=>(f.* ms.^2 .+ 0.5*v), 1=>fill(-0.25v,length(ms)-1), -1=>(fill(-0.25v,length(ms)-1))) .* csl
   tvals,tvecs = eigen(Matrix(ur(mc)*ht*ur(mc)))
   tvals = tvals[1:mm+1]
   tvecs = tvecs[:,1:mm+1]
   hr = hrotest([1.75;1.25;0.125;.08],n)
   rho = -2*f*csl*0.5.* kron((tvecs' * Diagonal(ms) * tvecs), Diagonal(collect(-n:n)))
   #@show rho
   out = dropzeros!(sparse(kron(I(mm+1),hr) + rho + Diagonal(kron(tvals,ones(2n+1)))))
   out = kron(I(mm+1),ur(n)) * out * kron(I(mm+1),ur(n))
   return dropzeros!(out) ./csl
end

function build_h2(n,mc,mm,σ)
   f = 5
   v = 200
   ms = collect(-(3mc+σ):3:(3mc+σ))
   ht = spdiagm(0=>(f.* ms.^2 .+ 0.5*v), 1=>fill(-0.25v,length(ms)-1), -1=>(fill(-0.25v,length(ms)-1))) .* csl
   tvals,tvecs = eigen(Matrix(ht))
   tvals = tvals[1:mm+1]
   tvecs = tvecs[:,1:mm+1]
   hr = hrotest([1.75;1.25;0.125;.08],n)
   rho = -2*f*csl*0.05.* kron((tvecs' * Diagonal(ms) * tvecs), Diagonal(collect(-n:n)))
   out = dropzeros!(sparse(kron(I(mm+1),hr) + rho + kron(Diagonal(tvals),I(2n+1))))
   out = kron(I(mm+1),ur(n)) * out * kron(I(mm+1),ur(n))
   return dropzeros!(out) ./csl
end

function build_h3(n,mc,mm,σ)
   f = 5
   v = 200
   ms = collect(-(3mc+σ):3:(3mc+σ))
   tbase = spdiagm(0=>(f.* ms.^2 .+ 0.5*v), 1=>fill(-0.25v,length(ms)-1), -1=>(fill(-0.25v,length(ms)-1))) .* csl
   for k in -n:n
   end
end

function twostg_op(prm::Float64,j::Real,s::Real,qns::Array{Int,2},ms::Array{Int},
                  a::Int,b::Int,c::Int,d::Int,e::Int,f::Int,g::Int,h::Int,
                  jp::Int,tvecs::Array{Float64,2})::SparseMatrixCSC{Float64, Int64}
   out = 0.25*prm*rsr_op(j,s,qns,a,b,c,d,e,f,jp)
   if (g≠0)||(h≠0)||(jp≠0)
      out = kron(tvecs' * tor_op(ms,g,h,jp) * tvecs,out)
   else
      out = kron(I(size(tvecs,2)), out)
   end
   return dropzeros!(tplus!(out))
end
function twostg_op(prm::Float64,j::Real,s::Real,qns::Array{Int,2},
            ms::Array{Int},plist::Array{Int},tvecs::Array{Float64,2})::SparseMatrixCSC{Float64, Int64}
   #1/2 from a+a', 1/2 from np^0 + nm^0
   out = twostg_op(prm,j,s,qns,ms,plist[1],plist[2],plist[3],plist[4],plist[5],
                     plist[6],plist[7],plist[8],plist[9],tvecs)
   return out
end

function hjstage(sof,cdf::Array,cdo::Array,j,s,nf,tvals,tvecs,ms)::SparseMatrixCSC{Float64, Int64}
   qns = qngen(j,s)
   #@show j
   #ms = msgen(nf,mc,σ)
   ℋ = hrot2(sof[1:4],qns) 
   if s==0.5
      ℋ .+= hsr(sof[5:9],j,s,qns)
   elseif s≥1
      ℋ .+= hsr(sof[5:9],j,s,qns)
      ℋ .+= hqu(sof[10:12],j,s,qns)
   end
   ℋ = kron(I(length(tvals)), ℋ)
   ℋ .+= kron(Diagonal(tvals), I(size(qns,1)))
   if (s≥0.5)&&(nf≠0)
   ℋ .+= kron(tvecs' *pa_op(ms,1)* tvecs, sof[14]*nz_op(qns,1) + sof[15]*npm_op(qns,1) + 
               sof[17]*sz_op(j,s,qns,1) + sof[18]*spm_op(j,s,qns,1))
   elseif (s==zero(s))&&(nf≠0)
      tpart = droptol!(sparse(tvecs' *pa_op(ms,1)* tvecs),1e-12)
      ℋ += kron(tpart, sof[14]*nz_op(qns,1) + sof[15]*npm_op(qns,1))
   else
   end
   for i in 1:length(cdf)
      ℋ .+= twostg_op(cdf[i],j,s,qns,ms,cdo[:,i], tvecs)
   end
   return dropzeros!(ℋ)
end

function twostg_diag(ctrl,sof,cdf,cdo,tvals,tvecs,ms,nf,mmax,mcalc,j,s,σ,vtm)
   #println("ℋ time for J=$j")
   H = hjstage(sof,cdf,cdo,j,s,nf,tvals,tvecs,ms)
   if true ∈ isnan.(H)
      @warn "FUCK!!! j=$j, σ=$σ, NaN in H"
   end
   U = kron(I(length(tvals)), ur(j,s))
   H = (U*H*U)
   vals, vecs = eigen!(Symmetric(Matrix(H)))

   perm = twostg_assign(vecs,j,s,mmax,vtm)
   vals = vals[perm]
   vecs = vecs[:,perm]
   
   vecs = U*vecs
   return vals, vecs
end

function twostg_calc(ctrl,prm,stg,cdo,nf,vtm,mcalc,mmax,jlist,s,sd,σ)
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   tormat, ms = htor2v2(sof[13:16],nf,mcalc,σ)
   tvals,tvecs = eigen(Symmetric(Matrix(tormat)))
   tvals = tvals[1:mmax+1]
   tvecs = tvecs[:,1:mmax+1]
   msd = sd*Int(mmax+1)
   vtd = Int(vtm+1)
   jmin = 0.5*iseven(sd)
   jmax = jlist[end]
   jfd = sd*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   outvals = zeros(Float64,jfd*vtd)
   outpasz = zeros(Float64,jfd*vtd)
   outquns = zeros(Int,jfd*vtd,6)
   outvecs = zeros(Float64,Int((2*jmax+1)*msd),jfd*vtd)
   for j in jlist #thread removed for troubleshooting purposes
#   @threads for j in jlist
      jd = Int(2.0*j) + 1
      sind, find = jvdest(j,s,vtm) 
      vals, vecs = twostg_diag(ctrl,sof,cdf,cdo,tvals,tvecs,ms,nf,mmax,mcalc,j,s,σ,vtm)
      #if j==1.0
      #   @show vals./csl
      #end
      outvals[sind:find] = vals
      outquns[sind:find,:] = qnlabv(j,s,nf,vtm,σ)
      outvecs[1:jd*msd,sind:find] = vecs
   end
   return outvals, outvecs, tvecs, outquns
end
function twostg_calc2(prm,stg,cdo,nf,ctrl,jlist)
#function tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
   s = ctrl["S"]
   mmax = ctrl["mmax"]
   mcalc = ctrl["mcalc"]
   vtm = ctrl["vtmax"]
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   vtd = Int(vtm+1)
   msd = Int((mmax+1)*(2s+1))
   jmin = 0.5*iseven((2*s+1))
   jmax = jlist[end,1]*0.5
   jfd = Int((2s+1)*sum(2. .*collect(Float64,jmin:jmax).+1.))
   σcnt = σcount(nf)
   fvls = zeros(Float64,jfd*vtd,σcnt)
   q = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):jmax) .+1))
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(mmax+1)),jfd*vtd,σcnt)
   tvecs = zeros(Float64,Int(2*mcalc+1),Int(mmax+1),σcnt)
   @time for sc in 1:σcnt
      σ = sc - 1
      tormat, ms = htor2v2(sof[13:16],nf,mcalc,σ)
      σvals,σvecs = eigen(Symmetric(Matrix(tormat)))
      σvals = σvals[1:mmax+1]
      σvecs = σvecs[:,1:mmax+1]
      tvecs[:,:,sc] = σvecs
      jmsd = Int(msd*(2*jmax+1))
      jsvd = Int(jfd*vtd)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      for j in jsublist
         jd = Int(2.0*j) + 1
         sind, find = jvdest(j,s,vtm) 
         fvls[sind:find,sc], fvcs[1:jd*msd,sind:find,sc] = twostg_diag(ctrl,
                           sof,cdf,cdo,σvals,σvecs,ms,nf,mmax,mcalc,j,s,σ,vtm)
         fqns[sind:find,:,sc] = qnlabv(j,s,nf,vtm,σ)
      end
   end
   return fvls, fvcs, tvecs, fqns
end

function vtfinder(svcs,jsd::Int,md::Int,vtmax)
   ovrlp = zeros(md,size(svcs,2))
   @simd for i in 1:md 
      ovrlp[i,:] = sum(svcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   vind = zeros(Int,size(svcs,1))
   cap = min(vtmax+4,md)
   for vi in 1:(vtmax+1)#cap #THIS HAS TO BE SERIAL DON'T SIMD THIS ONE FUTURE WES
      perm = sort(sortperm(ovrlp[vi,:], rev=true)[1:jsd])
      ovrlp[:,perm] .= 0.0
      ovrlp[vi,:] .= 0.0
      vind[perm] .= vi
   end
   return vind
end
function nfinderv4(svcs,vind,md,vtmax,jd,sd,ns,ni)
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
      nd = (2*ns[i]+1)*min((vtmax+1),md)
      #count = min(nd*(vtmax+1),nd*(md))
      #vlimit = min(vtmax+3,md) 
      #for v in 0:vlimit
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:nd])#[1:nd]
      nind[list[perm]] .= i
      ovrlp[:,perm] .= 0.0 
      #end
   end
   #println(nind)
   return nind
end
function twostg_assign(vecs,j,s,mmax,vtmax)
   md = mmax + 1
   jd = Int(2.0*j) + 1
   sd = Int(2.0*s) + 1
   ns, nd, ni, jsd = srprep(j,s)
   count = min(vtmax+4,md)
   svcs = abs.(vecs[:,1:jsd*count]).^2
   #determine vt state
   vind = vtfinder(svcs,jsd,md,vtmax)
   #determine N state
   nind = nfinderv4(svcs,vind,md,vtmax,jd,sd,ns,ni)
   #energy sort for K (done intrisnically)
   #permute to m style ordering (you didn't implement this....)
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for ng in 1:length(ns)
      nfilter = (nind .== ng)
      for v in 0:vtmax
         mg = mmax + vt2m(v) + 1
         frst = jsd*(mg-1) + ni[ng,1]
         last = jsd*(mg-1) + ni[ng,2]
         part = col[nfilter .* (vind .== (v+1))]
         part = part[keperm(ns[ng])]
         perm[frst:last] = part#col[filter][keperm(ns[ng])]
      end
   end
   perm = perm[perm .!= 0]
   return perm
end