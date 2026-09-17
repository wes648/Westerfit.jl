
eh(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
fhv(x::Float64,y::Int)::Float64 = □rt(x - eh(y))
dbfact(k::Int)::Int = prod(2k-1 : -2 : 1)
function fox(j::Number,k::Int)::Float64
   out = prod(doubled(j)-k : doubled(j)+k+1)*0.5^k
   out *= factorial(k)/dfact(k)
   return □rt(out)
end


################################ ROTATION |N,K⟩ OPERATORS ################################
function N2(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(mapreduce(x->fill(eh(x)^p, 2x+1), append!, ψ.N))
end

function Nz(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return spdiagm(mapreduce(x->(-x:x).^p, append!, ψ.N))
end

function Np(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64}
   if p ≠ 0 &&p ≤ ψ.lng
      ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])[1+p:end]
      ks = mapreduce(x->-x:x, vcat, ψ.N)[1+p:end]
      out = ones(length(ks))
      out = prod(fh.(ns,ks .- collect(1:p)'),dims=2)[:]
      out = spdiagm(-p=>out)
   elseif p ≠ 0 && p > ψ.lng
      out = spzeros(ψ.lng,ψ.lng)
   else
      out = spdiagm(ones(ψ.lng))
   end
   return out
end
function Npm(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
	return tplus!(Np(ψ,p))
end

Nx(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64} = tplus!(0.5.*Np(ψ,1))^p
function Ny2(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
   out = 0.5 .* (N2(ψ,1,0) - Nz(ψ,2,0))
   out -= 0.25 .* Npm(ψ,2,0)
   return out^p
end

function TN(ψ::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int64}
   # (-)^(N'-K') * wig3j(N',k,N, -K',q,K) * fox(N)
   out = spzeros(ψ.lng,ψ.lng)
   nds = nindsgen(ψ.N)
   for i ∈ 1:length(ψ.N)
      ks = -ψ.N[i]:ψ.N[i]
      blck = view(out,nds[i],nds[i])
      fact = fox(doubled(ψ.N[i]))
      blck[diagind(blck,q)] = fact .*
                              powneg1(ψ.N[i] .- ks[1:end-q]') .* 
                              wig3j(ψ.N[i],k,ψ.N[i], -ks[1:end-q]',q,ks[1:end-q])
      if !iszero(q)
      blck[diagind(blck,-q)]= fact .*
                              powneg1(k+q+ψ.N[i] .- ks[1:end-q]') .* 
                              wig3j(ψ.N[i],k,ψ.N[i], -ks[1:end-q]',-q,ks[1:end-q])
      end # q if
   end # n loop
   return out
end

################################## SPIN |S,Σ⟩ OPERATORS ##################################
function S2(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
   return spdiagm( fill(eh(ψ.S)^p,ψ.lng) )
end

function ts_fact(sf::Float64,j::Float64,s::Float64,nb::Int,nk::Int,k::Int)::Float64
   sf*jnred(nb,nk)*wig6j(s,nb,j, nk,s,k)
end
function wigeck_elem(x::Int,nb::Int,nk::Int,k::Int,q::Int,fac::Float64)::Float64
   fac*wig3j(nb,k,nk,-x-q,q,x)*powneg1(x+q)
end
function TS(ψ::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   @assert abs(q) ≤ k "T^k_q(S) component must be ≤ rank!"
   out = spzeros(ψ.lng,ψ.lng)
   nds = nindsgen(ψ.N)
   sfact = fox(ψ.S,k)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+k,length(ψ.N))
      fac = ts_fact(sfact,ψ.J,ψ.k, ψ.N[j], ψ.N[i],k)
      p = -q - ψ.N[j] + ψ.N[i] # -q - Δ
      blck = view(view(out,nds[j],nds[i]), diagind(dgen(ψ.N[j]),dgen(ψ.N[i]), p))
      kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
      map!(x::Int-> wigeck_elem(x,ψ.N[j],ψ.N[i],2,q,fac), blck ,kl)
      if q>0
         p = q - ψ.N[j] + ψ.N[i]
         blck = view(view(out,nds[j],nds[i]), diagind(view(out,nds[j],nds[i]), p))
         kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
         map!(x::Int-> wigeck_elem(x,ψ.N[j],ψ.N[i],2,q,fac*powneg1(q)), blck ,kl)
      end # q if
   end # i,j
   return Symmetric(out, :L)
end

function Sz(ψ::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   dropzeros!(TS(ψ,1,0)^p)
end
function Spm(ψ::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   dropzeros!(TS(ψ,k,k) .* (powneg1(k)*√2^k))
end


function hq_elem(x::Int,nb::Int,nk::Int,q::Int,fac::Float64)::Float64
   fac*wig3j(nb,2,nk,-x-q,q,x)*powneg1(x+q)
end
hq_sfac(j::Float64,s::Float64)::Float64 = 0.25*powneg1(j + s + 1.0) / wig3j(s,2,s, -s,0,s)
function hq_bfac(j::Float64,s::Float64,nb::Int,nk::Int,sf::Float64)::Float64
   wig6j(j, s,  nb,
         2,  nk, s)*jnred(nb,nk)*sf*powneg1(nb+nk)
end
function T2Q(ψ::RPsi,q::Int,p::Int)::SparseMatrixCSC{Float64,Int}
   @assert abs(q)≤2 "(T²(Q)⋅T²(V))_q only supports |q| ≤ 2"
   @assert ψ.S ≥ 1 "Spin must be at least 1 for quadrupole terms!"
   out = spzeros(ψ.lng,ψ.lng)
   nds = nindsgen(ψ.N)
   sfact = hq_sfac(ψ.J,ψ.S)
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+2,length(ψ.N))
      fac = hq_bfac(ψ.J,ψ.S, ψ.N[j],ψ.N[i], sfact)
      p = -q - ψ.N[j] + ψ.N[i] # -q - Δ
      blck = view(view(out,nds[j],nds[i]), diagind(dgen(ψ.N[j]),dgen(ψ.N[i]), p))
      kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
      map!(x::Int-> wigeck_elem(x,ψ.N[j],ψ.N[i],2,q,fac), blck ,kl)
      if q>0
         p = q - ψ.N[j] + ψ.N[i]
         blck = view(view(out,nds[j],nds[i]), diagind(view(out,nds[j],nds[i]), p))
         kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
         map!(x::Int-> wigeck_elem(x,ψ.N[j],ψ.N[i],2,-q,fac*powneg1(q)), blck ,kl)
      end # q if
   end # i,j
   return Symmetric(out, :L)
end

########################### SPIN-ROTATION |J,S,N,K⟩ OPERATORS ############################
function NS(ψ::RPsi,p::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   out = spdiagm(mapreduce(x-> fill( (0.5*(eh(ψ.J)-eh(ψ.S)-eh(x)))^p, 2x+1), vcat,ψ.N))
   return out
end

function srelem(x::Int,pr::Float64,nb::Int,nk::Int,l::Int,q::Int)::Float64
   pr*wig3j(nb,l,nk,-x-q,q,x)*powneg1(x)
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
function sr_blck_fac(j::Float64,s::Float64,nb::Int,nk::Int,sf::Float64)::Float64
   jsred(j,s,nb,nk)*nsred2(nb,nk)*sf
end

function TNS(ψ::RPsi, k::Int, q::Int)::SparseMatrixCSC{Float64,Int}
   @assert k ≤ 2 "T^k_q (N,S) only implemented up to rank-2"
   @assert abs(q) ≤ k "T^k_q (N,S) |q| cannot exceed k"
   @assert !isone(k) "Julia hasn't implemented T^1_q (N,S). pls email her if you need it"
if iszero(k)
   out = NS(ψ,1,0) .* -3^-0.5
else   
   out = spzeros(ψ.lng,ψ.lng)
   nds = nindsgen(ψ.N)
   sfact = √3*nred(ψ.S)*powneg1(ψ.J + ψ.S)
#   for i ∈ 1:length(ψ.N), j ∈ max(1,i-1):min(i+1,length(ψ.N))
   for i ∈ 1:length(ψ.N), j ∈ i:min(i+1,length(ψ.N))
#      nb = ψ.N[j]; nk = ψ.N[i]#; Δ = nb - nk
      fac = sr_blck_fac(ψ.J,ψ.S,ψ.N[j],ψ.N[i],sfact)
      p = -q - ψ.N[j] + ψ.N[i] # -q - Δ
      blck = view(view(out,nds[j],nds[i]), diagind(dgen(ψ.N[j]),dgen(ψ.N[i]), p))
      kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
      map!(x::Int-> srelem(x,fac,ψ.N[j],ψ.N[i],k,q), blck ,kl)
      if q>0
         p = q - ψ.N[j] + ψ.N[i] # -q - Δ
         blck = view(view(out,nds[j],nds[i]), diagind(dgen(ψ.N[j]),dgen(ψ.N[i]), p))
         kl = (-ψ.N[i]:ψ.N[i])[(1:length(blck)).+( p>0 ? p : 0)]
         map!(x::Int-> srelem(x,fac*powneg1(q),ψ.N[j],ψ.N[i],k,-q), blck ,kl)
      end # q
   end # i,j
   out = Symmetric(out, :L)
end # k if
   return out
end # func
################################# TOTAL |J,Ω⟩ OPERATORS ##################################

################################ TORSION |m,σ⟩ OPERATORS #################################
function Pt(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   if iszero(ψ.σ)
      out = map(x -> abs(x)^p, ψ.ms)
      if iseven(p)
         out = sparse(1:ψ.l, 1:ψ.l, out)
      else isodd(p)
         out = sparse(1:ψ.l, ψ.l:-1:1, out)
      end
   else
      out = sparse(1:ψ.l, 1:ψ.l, ψ.ms .^p)
   end
   return out
end

function cost(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   p = floor(Int, p/(ψ.nf * (1+iseven(ψ.nf)) ))
   out = spdiagm(p=>fill(0.5,ψ.l-p),-p=>fill(0.5,ψ.l-p))
   if iszero(ψ.σ)
      u = ul(ψ.l)
      out = dropzeros!(sand(out,u))
   end
   #torsetter!(ψ,tid,out)
   return out
end
function vnct(ψ::TPsi,p::Int,tid::Int)::SparseMatrixCSC{Float64, Int}
   p = floor(Int, p/(ψ.nf * (1+iseven(ψ.nf)) ))
   out = spdiagm(0=>fill(0.5,ψ.l),p=>fill(-0.25,ψ.l-p),-p=>fill(-0.25,ψ.l-p))
   if iszero(ψ.σ)
      u = ul(ψ.l)
      out = dropzeros!(sand(out,u))
   end
   #torsetter!(ψ,tid,out)
   return out
end

function sint(ψ::TTPsi,p::Int,tid::Int)::SparseMatrixCSC{ComplexF64, Int}
   l = 2ψ.mc+1
   out = spdiagm(p=>fill(0.5im,l-p),-p=>fill(-0.5im,l-p))
   if iszero(ψ.σs[tid])
      u = ur(ψ.mc)
      out = dropzeros!(u*out*u)
   end
   torsetter!(ψ,tid,out)
   return out
end

Pα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,1)
Pβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,2)
Pγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = Pt(ψ,p,3)

cosα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cost(ψ,p,1)
cosβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cost(ψ,p,2)
cosγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = cost(ψ,p,3)
vncα(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnct(ψ,p,1)
vncβ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnct(ψ,p,2)
vncγ(ψ::TPsi,p::Int,q::Int)::SparseMatrixCSC{Float64, Int} = vnct(ψ,p,3)

sinα(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sint(ψ,p,1)
sinβ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sint(ψ,p,2)
sinγ(ψ::TTPsi,p::Int,q::Int)::SparseMatrixCSC{ComplexF64, Int} = sint(ψ,p,3)

# https://doi.org/10.1103/PhysRevA.80.042513
# hirota 2.5.34
# edmonds 3.7.8
function Tμ(ψb,ψk,l,q)

end
