relem(j,k) = (0.5^k)*√( prod( 2j-k+1 : 2j+k+1 ) * factorial(k) / df(2k-1) )

#fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))

hr2on(ns,ks,bk::Float64,bn::Float64) = @. bn*eh(ns) + bk*ks^2 
hr2of1(ns,ks,dab::Float64) = @. dab*(ks-0.5)*fh(ns,ks-1)
hr2of2(ns,ks,bpm::Float64) = @. bpm*fh(ns,ks-1)*fh(ns,ks-2)

hr2onv(ns,ks,bk::Float64,bn::Float64) = bn*ns + bk*ks^2 
hr2of1v(ns,ks,dab::Float64) = dab*(ks-0.5)*fhv(ns,ks-1)
hr2of2v(ns,ks,bpm::Float64) = bpm*fhv(ns,ks-1)*fhv(ns,ks-2)

function hrot2(pr::Vector{Float64},ψ::RPsi)::SparseMatrixCSC{Float64, Int64}
   ns = nsgen(ψ.N)
   ks = ksgen(ψ.K)
   out = spdiagm(hr2on(ns,ks,pr[1],pr[2]))
   out[diagind(out,-1)] .= hr2of1(ns[2:end],ks[2:end], pr[4])
   out[diagind(out,-2)] .= hr2of2(ns[3:end],ks[3:end], pr[3])
   return dropzeros!(out)
end
function hrot2v(pr::Vector{Float64},ψ::RPsi)::SparseMatrixCSC{Float64, Int64}
   n2 = n2gen(ψ.N)
   ks = reduce(vcat, ψ.K)
   out = spdiagm(hr2onv.(n2,ks,pr[1],pr[2]))
   out[diagind(out,-1)] .= hr2of1v.(n2[2:end],ks[2:end], pr[4])
   out[diagind(out,-2)] .= hr2of2v.(n2[3:end],ks[3:end], pr[3])
   return dropzeros!(out)
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
   out = wig3j.(nb,l,nk,-kl,q,kl.-q) 
   out .*= powneg1.(kl)
   out .*= pr
   return out
end

function hsr(pr::Array{Float64},ψ::RPsi)::SparseMatrixCSC{Float64,Int}
   #pr = [T0_0 T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = √3*nred(S)*powneg1(J+S)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+1,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; ks = ψ.K[i]; Δ = nb - nk
      blck = view(out,nds[j],nds[i])
      frac = jsred(J,S,nb,nk)*nsred2(nb,nk)*sfact
      for p ∈ (-2-Δ):Δ
         q = Δ+p
         dest = diagind(blck,p)
         kl = ks[(1:length(dest)).+δi(1,p)]
         #the q in the phase factor is for T2_1 = -T2_-1
         blck[dest] = srelem(pr[2+abs(q)]*frac*powneg1(δi(q,-1)),nb,nk,kl,2,q)
      end#p loop
   end
   dropzeros!(out)
   #out .*= nred(S)*powneg1(J+S)
   out[diagind(out)] .+= -√(1/3)*pr[1] .* ns_el3(J,S,ψ.N)
   return out
end

function qured(j,s,nb,nk)
   return jnred(nb,nk)*wig6j(j, s,nb,
                             2,nk, s)
end
function quelem(pr::Float64,nb::Int,nk::Int,
         kl::UnitRange{Int},q::Int)::Vector{Float64}
   out = wig3j.(nb,2,nk,-kl,q,kl.-q) 
   out .*= powneg1.(kl)
   out .*= pr
end

function hqu(pr::Array{Float64},ψ::RPsi)::SparseMatrixCSC{Float64,Int}
   #pr = [T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = 0.25*inv(wig3j(S,2,S, -S,0,S))*powneg1(J+S+1)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+2,length(ψ.N))
      nk = ψ.N[j]; nb = ψ.N[i]; ks = ψ.K[i]; Δ = nb - nk
      blck = view(out,nds[j],nds[i])
      fac = qured(J,S,nb,nk)*powneg1(nb+nk+1)*sfact
      #rng = 
   Threads.@threads for q ∈ -2: 2*(1-δi(Δ,0))
         p = Δ+q
         dest = diagind(blck,p)
         kl = ks[(1:length(dest)).+δi(1,p)]
         #the q in the phase factor is for T2_1 = -T2_-1
         blck[dest] .= quelem(pr[1+abs(q)]*fac*powneg1(δi(q,-1)),nb,nk,kl,q)
      end
   end
   dropzeros!(out)
   return out
end

function htor(pr::Array{Float64},mc::Int,ψ::TPsi)::SparseMatrixCSC{Float64,Int}
   out = spzeros(ψ.lng,ψ.lng)
   for i ∈ 1:length(ψ.nf)
      r = length(ψ.nf) - i
      out += kron(I((2mc+1)^r), 
                 htorhc(pr[9+4i],pr[12+4i], mc, ψ.ms[i], ψ,σ[i]), 
                 I((2mc+1)^(i-1)) )
   end
   return out
end

function htorhc(f,v,mc,ms,σ)
   out = sparse(1:2mc+1, 1:2mc+1, map(x -> f*x^2 + v, ms))
   if isodd(nf)&&iszero(σ)
      out[diagind(out,-1)[1:mc-1]] .= -0.5v
      out[diagind(out,-1)[mc]] = -√0.5 * v
      out[diagind(out,-1)[mc+1:end]] .= -0.5v   
   elseif iseven(nf)&&iszero(σ)
      out[mc,mc] += 0.5v
      out[mc+2,mc+2] -= 0.5v
      out[diagind(out,-2)[1:mc-2]] .= -0.5v
      out[diagind(out,-2)[mc]] = -√0.5 * v
      out[diagind(out,-2)[mc+2:end]] .= -0.5v   
   else
      out[diagind(out, 1+iseven(nf))] .= -0.5*v
      out[diagind(out,-1-iseven(nf))] .= -0.5*v      
   end
   return out
end

function rhorho(rz,rx,ns,ks)::SparseMatrixCSC{Float64,Int}
   #shift is 4*topid
   out = spdiagm(rz .* ks)
   out[diagind(out,-1)] .= rx .* fh.(ns,ks[1:end-1])#<ψ|Nxψ> has an intrinsict factor of 1/2, cancels with -2Fρx
end

function rhomat(ctrl,pr::Array{Float64},ψ::Psi;tvecs=[1.0],otvcs=[1.0])::SparseMatrixCSC{Float64,Int}
   ns = nsgen(ψ.N)
   ks = ksgen(ψ.K)
   out = rhorho(pr[13]*pr[14],pr[13]*pr[15],ns,ks)
   stgchk = ctrl.stages > 2
   for i ∈ 2:length(ψ.T.nf)
      tpart = pa_op(ψ.T,2)
      if stgchk
         tpart = sand(tpart,sparse(otvcs[:,:,i])) #' * tpart * otvcs[:,:,i]
      end
      part = kron(tpart, rhorho(pr[10+4i],pr[11+4i],ns,ks))
      out = kron(sparse(I,size(part,1),size(part,1)), out) + kron(part, sparse(I,size(out,1),size(out,1)))
   end
   if ctrl.stages > 1
      out = sand(out,tvecs)
   end
   return out#sparse(out)
end

function htsr2_1stg(ctrl::Controls,pr::Array{Float64},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   H = hrot2(prm[1:4],ψ)
   if ctrl.S ≥1.0
      H += hsr(prm[5:8],ψ.R.J,ψ.R.S,ψ.R) + hqu(prm[9:11],ψ.R.J,ψ.R.S,ψ.R)
   elseif ctrl.S ==0.5
      H += hsr(prm[5:8],ψ.R.J,ψ.R.S,ψ.R)
   end
   H = kron(I(ψ.T.lng),H)
   H += kron( htor(pr, ctrl.mc, ψ.T), I(ψ.R.lng) )
   H += rhomat(ctrl,vals,ψ)
   return Symmetric(H,:L)
end
function htsr2_2stg(ctrl::Controls,pr::Array{Float64},ψ::Psi,tvals::Array{Float64},tvecs::Array{Float64,2})::SparseMatrixCSC{Float64,Int}
   H = hrot2(prm[1:4],ψ)
   if ctrl.S ≥1.0
      H += hsr(prm[5:8],ψ.R.J,ψ.R.S,ψ.R) + hqu(prm[9:11],ψ.R.J,ψ.R.S,ψ.R)
   elseif ctrl.S ==0.5
      H += hsr(prm[5:8],ψ.R.J,ψ.R.S,ψ.R)
   end
   H = kron(I(ψ.T.lng),H)
   H += kron( Diagonal(tvals), sparse(I,ψ.R.lng,ψ.R.lng) )
   H += rhomat(ctrl,vals,ψ,tvecs=tvecs)
   return Symmetric(H,:L)
end


#= This is the graveyard

function hsr_coarse(pr::Array{Float64},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])
   ks = reduce(vcat, ψ.K)
   out  = @. wig3j(ns', 2, ns, -ks', -2, ks)
   out += @. wig3j(ns', 2, ns, -ks', -1, ks)
   out += @. wig3j(ns', 2, ns, -ks',  0, ks)
   out += @. wig3j(ns', 2, ns, -ks',  1, ks)
   out += @. wig3j(ns', 2, ns, -ks',  2, ks)
   out[abs.(ns' .- ns) .> 1] .*= 0.0
   return out
end
function sz_coarse(pr::Array{Float64},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])
   ks = reduce(vcat, ψ.K)
   out = @. wig3j(ns', 2, ns, -ks',  0, ks)
   out[abs.(ns' .- ns) .> 1] .*= 0.0
   return out
end
function sp_coarse(ψ::Psi)::SparseMatrixCSC{Float64,Int}
   ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])
   ks = reduce(vcat, ψ.K)
   out = @. wig3j(ns', 1, ns, -ks',  1, ks)
   out[abs.(ns' .- ns) .> 1] .*= 0.0
   return out
end

function hqu_coarse(pr::Array{Float64},ψ::Psi)::SparseMatrixCSC{Float64,Int}
   ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])
   ks = reduce(vcat, ψ.K)
   out  = @. wig3j(ns', 2, ns, -ks', -2, ks)
   out += @. wig3j(ns', 2, ns, -ks', -1, ks)
   out += @. wig3j(ns', 2, ns, -ks',  0, ks)
   out += @. wig3j(ns', 2, ns, -ks',  1, ks)
   out += @. wig3j(ns', 2, ns, -ks',  2, ks)
   out[abs.(ns' .- ns) .> 2] .*= 0.0
   return out
end


function hrot2_old(pr::Vector{Float64},ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   out = spzeros(size(ψ.N,1),size(ψ.N,1))
   #p0 = hr2on(ns,ks,pr[1],pr[2])
   out[diagind(out)] .= hr2on(ψ.N,ψ.K,pr[1],pr[2])
   #p1 = hr2of1(ns[2:end],ks[2:end], pr[4])
   out[diagind(out,1)] .= hr2of1(ψ.N[2:end],ψ.K[2:end], pr[4])
   #p2 = hr2of2(ns[3:end],ks[3:end], pr[3])
   out[diagind(out,2)] .= hr2of2(ψ.N[3:end],ψ.K[3:end], pr[3])
   #out = spdiagm(0=>p0,1=>p1,2=>p2)
   return dropzeros!(out)
end

function nsred(l::Int,nb::Int,nk::Int)::Float64
   if (abs(l)==1)&&(nb==nk)
      out = 0.0
   else
   out = 0.5*( 
   √(nk*(nk+1)*(2*nk+1))*
      wig6j( 1, 1, l,
            nb,nk,nk) + #*powneg1(l) + 
   √(nb*(nb+1)*(2*nb+1))*
      wig6j( 1, 1, l,
            nk,nb,nb))
   end
   return out
end

function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)::Float64
   return pr*wig3j( nb,l,nk,
                   -kb,q,kk)*√(2*l+1)*
       nsred(l,nb,nk)*jsred(j,s,nb,nk)*powneg1(nb-nk-kb)
end
function hsr(pr::Array{Float64},j,s,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   le = length(ψ.K)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   ts = SA[0  1 1 2  2 2  2 2;
           0 -1 1 0 -1 1 -2 2;
           1  2 2 3  4 4  5 5;
           1  1 1 1 -1 1  1 1]
   for i in 1:8
      tv = ts[:,i]
      prm = pr[tv[3]]*tv[4]
      if prm ≠ 0.0
      for a in 1:le, b in a:le
         nb = ψ.N[b]
         kb = ψ.K[b]
         nk = ψ.N[a]
         kk = ψ.K[a]
         if abs(nb-nk)≤1 && (tv[2]+kk-kb)==0
            out[b,a] += srelem(prm,tv[1],tv[2], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#sr ind for loop
      end#prm chck if
   end#sr term for loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(j+s)
   return out
end
function wiginv(s::Real)::Float64
   return s<one(s) ? 0.0 : inv(wig3j( s,2,s,
                                     -s,0,s))
end
function qulm(pr,q,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   return pr*qured(j,s,nb,nk)*
#             δ(nb,nk)* #This line can be used to emulate the perturbative 
             wig3j( nb, 2,nk,
                   -kb, q,kk)*powneg1(nb+nk-kb)
end
function hqu(pr,j,s,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   le = size(ψ.K,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   tq = SA[0 -1 1 -2 2; 
           1  2 2  3 3;
           1 -1 1  1 1]
   for i in 1:5
      tv = tq[:,i]
      prm = pr[tv[2]]*tv[3]
      if prm ≠ 0.0
      for a in 1:le, b in a:le
         nb = ψ.N[b]
         kb = ψ.K[b]
         nk = ψ.N[a]
         kk = ψ.K[a]
         if abs(nb-nk)≤2 && (tv[1]+kk-kb)==0
            out[b,a] += qulm(prm,tv[1], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#qu ind for loop
      end#prm chck if
   end#qu term for loop
   dropzeros!(out)
   out .*= 0.25*wiginv(s)*powneg1(j+s)
   #@show out
   return out
end

=#
