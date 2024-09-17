"""
This file is for testing variations on the theme of multiple torsional stages
"""

include("/home/wes/files/westerfit/src/main.jl")

function build_h1(n,mc,σ)
   f = 5
   v = 200
   ms = collect(-(3mc+σ):3:(3mc+σ))
   ht = spdiagm(0=>(f.* ms.^2 .+ 0.5*v), 1=>fill(-0.25v,length(ms)-1), -1=>(fill(-0.25v,length(ms)-1))) .* csl
   hr = hrotest([1.75;1.25;0.125;.08],n)
   rho = -2*f*csl*0.05 .* kron(Diagonal(ms),Diagonal(collect(-n:n)))
   out = dropzeros(kron(ht,I(2n+1)) + kron(I(length(ms)),hr) + rho)
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
   rho = -2*f*csl*0.05.* kron((tvecs' * Diagonal(ms) * tvecs), Diagonal(collect(-n:n)))
   #@show rho
   out = dropzeros(sparse(kron(I(mm+1),hr) + rho + Diagonal(kron(tvals,ones(2n+1)))))
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
   out = dropzeros(sparse(kron(I(mm+1),hr) + rho + kron(Diagonal(tvals),I(2n+1))))
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
      out = kron(I(size(ms,1)), out)
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
      ℋ += kron(tvecs' *pa_op(ms,1)* tvecs, sof[14]*nz_op(qns,1) + sof[15]*npm_op(qns,1))
   else
   end
   for i in 1:length(cdf)
      ℋ .+= twostg_op(cdf[i],j,s,qns,ms,cdo[:,i], tvecs)
   end
   return dropzeros!(ℋ)
end

function twostg_diag(ctrl,sof,cdf,cdo,tvals,tvecs,ms,nf,mcalc,j,s,σ,vtm)
   #println("ℋ time for J=$j")
   H = hjstage(sof,cdf,cdo,j,s,nf,tvals,tvecs,ms)
   if true ∈ isnan.(H)
      @warn "FUCK!!! j=$j, σ=$σ, NaN in H"
   end
   U = kron(I(length(tvals)), ur(j,s))
   H = (U*H*U)
   vals, vecs = eigen!(Symmetric(Matrix(H)))

   perm = twostg_assign(vecs,j,s,mcalc,vtm)
   vals = vals[perm]
   vecs = vecs[:,perm]
   
   vecs = U*vecs
   return vals, vecs
end

function twostg_calc(ctrl,prm,stg,cdo,nf,vtm,mcalc,mlim,jlist,s,sd,σ)
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   tormat, ms = htor2v2(sof[13:16],nf,mcalc,σ)
   tvals,tvecs = eigen(Symmetric(Matrix(tormat)))
   tvals = tvals[1:mlim+1]
   tvecs = tves[:,1:mlim+1]
   msd = sd*Int(mlim+1)
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
      vals, vecs = twostg_diag(ctrl,sof,cdf,cdo,tvals,tvecs,ms,nf,mcalc,j,s,σ,vtm)
      outvals[sind:find] = vals
      outquns[sind:find,:] = qnlabv(j,s,nf,vtm,σ)
      outvecs[1:jd*msd,sind:find] = vecs
   end
   return outvals, outvecs, outquns
end
