#This is the hard coded hamiltonian file!
# it is intended to contain better optimized versions of the commonly called operations

function prm_proc(prm::Vector{Float64},lnf)::Vector{Float64}
   out = zero(prm)
   out[2] = 0.5*(prm[2] + prm[3])
   out[1] = prm[1] - out[2]
   out[3] = 0.25*(prm[2] - prm[3])
   out[1:3] .*= 0.5
   out[4:end] = prm[4:end]
   for i ∈ 1:lnf
      out[hccount + 4i - 3] *= csl
      out[hccount + 4i] *= csl*0.5
   end
   return out
end

hr2onv(ns,ks,bk::Float64,bn::Float64) = bn*ns + bk*ks^2 
hr2of1v(ns,ks,dab::Float64) = dab*(ks-0.5)*fhv(ns,ks-1)
hr2of2v(ns,ks,bpm::Float64) = bpm*fhv(ns,ks-1)*fhv(ns,ks-2)

n2gen(x::Int)::Vector{Float64} = fill(eh(x), 2x+1)

function hrot2_hc(prm,ns::UnitRange)::SparseMatrixCSC{Float64,Int}
   nv = mapreduce(n2gen, append!, ns)
   kv = mapreduce(x -> -x:x, vcat, ns)
   out = spzeros(size(nv,1),size(nv,1))
   out[diagind(out   )] .= hr2onv.(nv,kv,prm[1],prm[2])
   out[diagind(out,-1)] .= hr2of1v.(nv,kv, prm[4])[2:end]
   out[diagind(out,-2)] .= hr2of2v.(nv,kv, prm[3])[3:end]
   out[diagind(out, 1)] .= out[diagind(out,-1)]
   out[diagind(out, 2)] .= out[diagind(out,-2)]
   return dropzeros!(out)
end
function hrot2_hc!(out::SparseMatrixCSC{Float64,Int},prm::Vector{Float64},ns::UnitRange)#::SparseMatrixCSC{Float64,Int}
   nv = mapreduce(n2gen, append!, ns)
   kv = mapreduce(x -> -x:x, vcat, ns)
   out[diagind(out   )] .= hr2onv.(nv,kv,prm[1],prm[2])
   out[diagind(out,-1)] .= hr2of1v.(nv[2:end],kv[2:end], prm[4])
   out[diagind(out,-2)] .= hr2of2v.(nv[3:end],kv[3:end], prm[3])
   return nothing #dropzeros!(out)
end


function nindsgen(ns::UnitRange{Int})::Vector{UnitRange{Int}}
   ni = Vector{UnitRange{Int}}(undef,length(ns))
   ni[1] = 1:2*ns[1]+1
   for i in 2:length(ns)
      ni[i] = (ni[i-1][end]+1):(ni[i-1][end]+2ns[i]+1)
   end
   return ni
end
function nsred2(nb::Int,nk::Int)::Float64
   if nb==nk
   out = nred(nk)*wig6j(1,1,2,nk,nk,nk)
   else
   out = 0.5*(nred(nk)*wig6j(1,1,2, nb,nk,nk) +
      nred(nb)*wig6j(1,1,2, nk,nb,nb))
   end
   return out
end
jnred(j::Real,n::Real)::Float64 = √((2*j+1)*(2*n+1))
nred(n::Real)::Float64 = √(n*(n+1)*(2*n+1))
function jsred(j,s,nb::Int,nk::Int)::Float64
   return wig6j(nk, s, j,
                 s,nb, 1)*jnred(nb,nk)
end

function srelem(pr::Float64,nb::Int,nk::Int,
         kl::UnitRange{Int},l::Int,q::Int)::Vector{Float64}
   #out = wig3j.(nb,l,nk,-kl,q,kl.-q) 
   #out .*= powneg1.(kl)
   #out .*= pr
   out = map(x-> pr*wig3j(nb,l,nk,-x,q,x-q)*powneg1(x), kl)
   return out
end
function ns_el3(j,s,ns::UnitRange)::Vector{Float64}
   reduce(vcat, [fill(0.5*(eh(j) - eh(n) - eh(s)), 2n+1) for n ∈ ns ])
end

ns_el(p,j,s,ns) = mapreduce(n->fill(0.5*p*(eh(j)-eh(s)-eh(n)), 2n+1), append!, ns)
function hsr!(out::SparseMatrixCSC{Float64,Int},pr::Array{Float64},ψ::RPsi)::SparseMatrixCSC{Float64,Int}
   #pr = [T0_0 T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   #out = spzeros(ψ.lng,ψ.lng)
   sfact = √3*nred(S)*powneg1(J+S)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+1,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; Δ = nb - nk
      blck = view(out,nds[j],nds[i])
      frac = jsred(J,S,nb,nk)*nsred2(nb,nk)*sfact
      for p ∈ (-2-Δ):Δ
         q = Δ+p
         dest = diagind(blck,p)
         kl = (-nk:nk)[(1:length(dest)).+δi(1,p)]
         #the q in the phase factor is for T2_1 = -T2_-1
         blck[dest] .+= srelem(pr[2+abs(q)]*frac*powneg1(δi(q,-1)),nb,nk,kl,2,q)
      end#p loop
   end
   dropzeros!(out)
   out[diagind(out)] .+= ns_el(-√(1/3)*pr[1], J,S,ψ.N)
   return nothing
end

function qured(j,s,nb,nk)::Float64
   return jnred(nb,nk)*wig6j(j, s,nb,
                             2,nk, s)
end
function quelem(pr::Float64,nb::Int,nk::Int,
         kl::UnitRange{Int},q::Int)::Vector{Float64}
   #out = wig3j.(nb,2,nk,-kl,q,kl.-q) 
   #out .*= powneg1.(kl)
   #out .*= pr
   return map(x-> pr*wig3j(nb,2,nk,-x,q,x-q)*powneg1(x), kl)
end
function hqu!(out::SparseMatrixCSC{Float64,Int},pr::Array{Float64},ψ::RPsi)::SparseMatrixCSC{Float64,Int}
   #pr = [T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   sfact = 0.25*inv(wig3j(S,2,S, -S,0,S))*powneg1(J+S+1)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+2,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; Δ = nb - nk
      blck = view(out,nds[j],nds[i])
      fac = qured(J,S,nb,nk)*powneg1(nb+nk+1)*sfact
      for q ∈ -2: 2*(1-δi(Δ,0))
         p = Δ+q
         dest = diagind(blck,p)
         kl = (-nk:nk)[(1:length(dest)).+δi(1,p)]
         #the q in the phase factor is for T2_1 = -T2_-1
         blck[dest] .+= quelem(pr[1+abs(q)]*fac*powneg1(δi(q,-1)),nb,nk,kl,q)
      end
   end
   dropzeros!(out)
   return nothing #out
end

function htor2_hc(pr::Vector{Float64},ψ::TPsi)::SparseMatrixCSC{Float64,Int}
   if isodd(ψ.nf)
      temp = fill(-0.5*pr[4], ψ.l-1)
      out = spdiagm(0=>map(x-> pr[1]*x^2 + pr[4], ψ.ms),
         1=>temp, -1=>temp)
   else # iseven(ψ.nf)
      temp = fill(-0.5*pr[4], ψ.l-2)
      out = spdiagm(0=>map(x-> pr[1]*x^2 + pr[4], ψ.ms),
         2=>temp, -2=>temp)
   end
   if iszero(ψ.σ)
      out = sand(out, ul(ψ.l))
   end
   return out
end

function htsr2_hc(pr,wvs,ψ,σid)::SparseMatrixCSC{Float64,Int}
   # U'Pt U ⊗ (ρz Nz + ρx Nx + ηz Sz + ηx Sx)
   NzM = Diagonal{Float64}(mapreduce(x -> -x:x, vcat, ψ.R.N))
   NxM = Npm(ψ.R,1,0)
#   SzM =  Sz(ψ.R,1,0)
#   SxM =  Sx(ψ.R,1,0)
   if length(ψ.T.nfs)>1 && !isnothing(wvs.top)
      torpart = spzeros(ψ.T.l, ψ.T.l)
      for i ∈ 1:length(ψ.T.nfs)
         torpart = torsetter!(ψ.T,i, sand(Pt(ψ.T.tps[i],1,0), wvs.top[i].vecs[:,:, σ2ind(ψ.T.σs[i], i, ψ.T.nfs[i])]))
      end
   elseif length(ψ.T.nfs)>1 && isnothing(wvs.top)
      torpart = spzeros(ψ.T.l, ψ.T.l)
      for i ∈ 1:length(ψ.Tnfs)
         torpart += torsetter!(ψ,i, Pt(ψ.top[i],1,0) )
      end
   elseif length(ψ.Tnfs)>1
      torpart = Pt(ψ.top[1],1,0)
   else
      torpart = spzeros(1,1)
   end
   if !isnothing(wvs.ttp)
      torpart = sand(torpart, wvs.ttp.vecs[:,:,σid])
   end
   out = kron( torpart, pr[hccount+2]*NzM + pr[hccount+3]*NxM) #+ pr[hccount+5]*SzM + pr[hccount+6]*SxM)
   for i ∈ 2:length(ψ.T.nfs)
      out += kron( torpart, pr[hccount+6i-4]*NzM + pr[hccount+6i-3]*NxM )
         #+ pr[hccount+6i-1]*SzM + pr[hccount+6i]*SxM)
   end
   return out
end


#=
reminder for the SR functions
julia> out = spzeros(12,12); @benchmark hsr!($out,$prm,$ψ) # current
BenchmarkTools.Trial: 7632 samples with 1 evaluation per sample.
 Range (min … max):  297.443 μs …   6.432 ms  ┊ GC (min … max): 0.00% … 83.93%
 Time  (median):     628.759 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   649.538 μs ± 171.527 μs  ┊ GC (mean ± σ):  0.52% ±  2.43%

                     ▂▃▄▄▅▇██▇▆▅▄▄▃▂▁▂▁▂▁▁▁                     ▂
  ▄▅▅▅▅▅▅▅▃▅▃▄▃▁▆▅▆████████████████████████▇▇█▆▇▆▇▇▆▆▅▆▅▅▄▅▅▅▅▄ █
  297 μs        Histogram: log(frequency) by time       1.06 ms <

 Memory estimate: 34.75 KiB, allocs estimate: 687.

julia> @benchmark hsr($prm,$ψ) # old version as reference
BenchmarkTools.Trial: 6498 samples with 1 evaluation per sample.
 Range (min … max):  321.240 μs …   9.727 ms  ┊ GC (min … max): 0.00% … 90.38%
 Time  (median):     710.734 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   758.674 μs ± 273.698 μs  ┊ GC (mean ± σ):  0.70% ±  2.93%

             ▃▇██▇▆▄▃▃                                           
  ▂▂▂▂▂▃▃▄▄▅███████████▇▇▅▅▄▄▃▃▃▃▃▃▂▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂ ▄
  321 μs           Histogram: frequency by time         1.81 ms <

 Memory estimate: 38.84 KiB, allocs estimate: 657.

julia> @benchmark hsr($prm, $J, $S, $qns) # westerfit v1 version
BenchmarkTools.Trial: 8593 samples with 1 evaluation per sample.
 Range (min … max):  315.418 μs … 184.546 ms  ┊ GC (min … max):  0.00% … 72.87%
 Time  (median):     352.117 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   573.483 μs ±   4.172 ms  ┊ GC (mean ± σ):  13.69% ±  1.93%

  ▆█▆▅▄▃▃▃▃▂▂▂▁▂▂▁▁▁▁            ▁▁▁▁▁  ▁▁                      ▂
  █████████████████████████▇███████████████▇▇▇█▇▆▇▆▆▆▅▅▅▅▆▄▅▅▅▅ █
  315 μs        Histogram: log(frequency) by time        1.2 ms <

 Memory estimate: 106.28 KiB, allocs estimate: 2972.
=#

