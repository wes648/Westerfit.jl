"""
Welcome to be base Operators used in westerfit! 
This file is divided into 3 sections and to add a new Operator, all 3 must be
  edited.
Section 1: Functions. These act on the wavefuntions to produce the matrices and
must use the following type annotations:
 op(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int64}
This should ensure consistency with the program.
Stage 0 operators (those used in the intensity calculations) have the form of:
 op(ψb::Psi,ψk::Psi,p::Int)::SparseMatrixCSC{Float64,Int64}
This is because we can mixing of J values in transitions but not energy levels

Section 2: Operators. These are the actual objects used by the code. Look at the
type file to see the structure but the general set up should be:
   Nam = Op(1.0,rp=[1],rf=[function])
where Nam is the name that will be read from the input file and function is the
new Operator function. We are defining it as being raised to the 1st power with
a default parameter value of 1.0

Section 3: Dictionary. This is the list of Operators read by the code in order
to interpret the input file. It is wrapped inside a function so that it doesn't 
stay in memory after being used in the input reading.
Each entry should just be "Nam"=>Nam
"""

################################################################################
#####                         Section 1: Functions                         #####
################################################################################
eh(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
fhv(x::Float64,y::Int)::Float64 = □rt(x - eh(y))
ns_el(j,s,p,n)::Float64 = (0.5*eh(j) - eh(n) - eh(s))^p

n2gen(ns::UnitRange{Int})::Vector{Float64} = reduce(vcat, [fill(eh(n),2n+1) for n ∈ ns])
nsgen(ns::UnitRange{Int})::Vector{Int} = reduce(vcat, [fill(n,2n+1) for n ∈ ns])
ksgen(ks::Vector{UnitRange{Int}})::Vector{Int} = reduce(vcat, ks)

#These functions are cooked into enact_init as they are purely diagonal
#they are N^2a (NS)^b S^c Nz^d in this order
function n2(ns::UnitRange{Int},p::Int)::Vector{Float64}
   reduce(vcat, [fill(eh(n)^p,2n+1) for n ∈ ns])
end
function n2(v::Float64,ns::UnitRange{Int},p::Int)::Vector{Float64}
   reduce(vcat, [fill(v*eh(n)^p,2n+1) for n ∈ ns])
end
function ns(j::Float64,s::Float64,ns::UnitRange,p::Int)::Vector{Float64}
   reduce(vcat, [fill((eh(j) - eh(n) - eh(s))^p, 2n+1) for n ∈ ns ])
end
function ns_el3(j,s,ns::UnitRange)::Vector{Float64}
   reduce(vcat, [fill(0.5*(eh(j) - eh(n) - eh(s)), 2n+1) for n ∈ ns ])
end
nz(ks::Vector{UnitRange{Int}},p::Int)::Vector{Float64} = reduce(vcat, ks).^p

function np_old(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int64}
   ns = ψ.N[1+p:end]
   part = ones(length(ns))
   if p ≤ ψ.lng && p ≠ 0
      ks = ψ.K[1+p:end]
      part = ones(length(ks))
      for o in 1:p
         part .*= fh.(ns,ks.-o)
      end
      out = spzeros(ψ.lng,ψ.lng)
      out[diagind(out,-p)] = part
   elseif p > ψ.lng && p ≠ 0
      out = spzeros(ψ.lng,ψ.lng)
      #out[diagind(out,-p)] = part
   else
      out[diagind(out)] .= 1.0 # = spdiagm(ones(ψ.lng))
   end
   return out
end

function np(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64}
   if p ≤ ψ.lng && p ≠ 0
      ns = reduce(vcat, [fill(n,2n+1) for n ∈ ψ.N])[1+p:end]
      ks = reduce(vcat, ψ.K)[1+p:end]
      out = ones(length(ks))
      out = prod(fh.(ns,ks .- collect(1:p)'),dims=2)[:]
      out = spdiagm(-p=>out)
   elseif p > ψ.lng && p ≠ 0
      out = spzeros(ψ.lng,ψ.lng)
   else
      out = spdiagm(ones(ψ.lng))
   end
   return out
end
nm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64} = permutedims(np(ψ,p))
npm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64} = tplus!(np(ψ,p))
nx(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int64} = tplus!(0.5*np(ψ,1))^p
function iny(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
   if p≠0
      out = np_op(ψ,1)
      out .-= permutedims(out)
      return dropzeros!(out)
   else
      return spdiagm(ones(ψ.lng))
   end
end

function s2(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
   return spdiagm( fill(eh(ψ.S)^p,ψ.lng) )
end

function skqpart(j,s,k,q,nb,kb,nk,kk)::Float64
   return wig3j(nb,k,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,k)*jnred(nb,nk)*powneg1(-kb)
end
function s0part(j,s,nb,kb,nk,kk)::Float64
   return skqpart(j,s,1,0,nb,kb,nk,kk)::Float64
end

#sz(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = sz_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,p)
sz(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = sz_op_new(ψ)^p
function sz_op(j::Float64,s::Float64,ns::UnitRange{Int},
       ks::Vector{UnitRange{Int}},lng::Int,p::Int)::SparseMatrixCSC{Float64,Int}
   nds = nindsgen(ns)
   out = spzeros(lng,lng)
   for i ∈ 1:length(ns); j ∈ max(i-1,1):i
      nb = ns[j]; nk = ns[i]
      blck = view(out,nds[j],nds[i])
      for k ∈ -2:2
      blck[diagind(blck,nb-nk+k)] = @. wig3j(nb,1,nk,-ks[j],k,ks[i])*powneg1(-kb)
      end
      blck .*= wig6j(s,nb,j,nk,s,1)*jnred(nb,nk)
   end#n loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(s+j+1)
   out = Symmetric(out,:L)^p
   return out
end
function sz_op_new(ψ::RPsi)::SparseMatrixCSC{Float64,Int}
   #pr = [T0_0 T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = nred(S)*powneg1(J+S)
   for i ∈ 1:length(ψ.N), j ∈ max(i-1,1):min(i+1,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; ks = ψ.K[i]; Δ = nk - nb
      blck = view(out,nds[j],nds[i])
      frac = wig6j(S,nb,J,nk,S,1)*jnred(nb,nk)*sfact
      #p = Δ
      dest = diagind(blck,Δ)
      kl = ks[(1:length(dest)).+δi(1,Δ)]
      blck[dest] = srelem(frac,nb,nk,kl,1,0)
   end
   dropzeros!(out)
   return out
end
#sq(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,p)
sq(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = sq_op_new(ψ,p)
function sq_op(j::Float64,s::Float64,ns::UnitRange{Int},
       ks::Vector{UnitRange{Int}},lng::Int,p::Int)::SparseMatrixCSC{Float64,Int}
   nds = nindsgen(ns)
   out = spzeros(lng,lng)
   for i ∈ 1:length(ns); j ∈ max(i-p,1):min(i+p,length(ns))
      nb = ns[j]; nk = ns[i]
      blck = view(out,nds[j],nds[i])
     blck[diagind(blck,nb-nk+p)] = @. wig3j(nb,p,nk,-ks[j],p,ks[i])*powneg1(-kb)
      blck .*= wig6j(s,nb,j,nk,s,p)*jnred(nb,nk)
   end#n loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(s+j+p)
   return out
end
function sq_op_new(ψ::RPsi,q::Int)::SparseMatrixCSC{Float64,Int}
   #pr = [T0_0 T2_0 T2_1 T2_2]
   J = ψ.J
   S = ψ.S
   nds = nindsgen(ψ.N)
   out = spzeros(ψ.lng,ψ.lng)
   sfact = nred(S)*powneg1(J+S+q)
   for i ∈ 1:length(ψ.N), j ∈ max(i-q,1):min(i+q,length(ψ.N))
      nb = ψ.N[j]; nk = ψ.N[i]; ks = ψ.K[i]; Δ = nk - nb
      blck = view(out,nds[j],nds[i])
      frac = wig6j(S,nb,J,nk,S,q)*jnred(nb,nk)*sfact
      p = Δ+q
      dest = diagind(blck,p)
      kl = ks[(1:length(dest))] .+ Δ
      blck[dest] = srelem(frac,nb,nk,kl,q,-q)
   end
   dropzeros!(out)
   return out
end


sp(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = p*sq(ψ,p)
sm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = p*sq(ψ,-p)
spm(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = tplus!(p*sq(ψ,p))
function sx(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
#   out  = sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng, 1)
#   out -= sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,-1)
   out  = sq_op_new(ψ, 1)
   out -= sq_op_new(ψ,-1)
   out *= -√.5
   out ^= p
   return out
end
function isy(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int}
   out  = sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng, 1)
   out -= sq_op(ψ.J,ψ.S,ψ.N,ψ.K,ψ.lng,-1)
   out *= -√.5
   out ^= p
   return out
end

E(ψ::RPsi,p::Int)::SparseMatrixCSC{Float64,Int} = sparse(1.0I,ψ.lng,ψ.lng)

function μ_gen(ψb::RPsi,ψk::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
#This matrix is transposed FUCK FUCK FUCK apr 29
   out = spzeros(ψb.lng,ψk.lng)
   nbinds = nindsgen(ψb.N)
   nkinds = nindsgen(ψk.N)
   for a ∈ 1:length(ψb.N), b ∈ 1:length(ψk.N)
   nbind = nbinds[a] 
   nkind = nkinds[b]
   nb = ψb.N[a]
   nk = ψk.N[b]
   if abs(nb-nk) ≤ k
   factr = wig6j(  nb,ψb.J,ψb.S, 
                 ψk.J,  nk,   k)*
           jnred(ψb.J,ψk.J)*jnred(nb,nk)*powneg1(nb+nb+ψk.S+ψk.J + k)
   for j ∈ 1:length(nkind), i ∈ 1:length(nbind)
      kb = -nb + i - 1
      kk = -nk + j - 1
      if iszero(kk+q-kb)
         out[nbind[i],nkind[j]] = wig3j(nb,k,nk,-kb,q,kk)*factr
      end;end #k loop;if
   end;end#n loop;if
   dropzeros!(out)
   return out
end#function
function μzf(ψb::RPsi,ψk::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   return μ_gen(ψb,ψk,1,0)^k
end
function μxf(ψb::RPsi,ψk::RPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   return dropzeros!(√0.5 .*(μ_gen(ψb,ψk,1,-1) - μ_gen(ψb,ψk,1,1)))^k
end
function iμyf(ψb::Psi,ψk::Psi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   return dropzeros!(-√0.5 .*(μ_gen(ψb,ψk,1,-1) + μ_gen(ψb,ψk,1,1)))^k
end
function μyf(ψb::Psi,ψk::Psi,k::Int,q::Int)::SparseMatrixCSC{ComplexF64,Int}
   return dropzeros!(im*√0.5 .*(μ_gen(ψb,ψk,1,-1) + μ_gen(ψb,ψk,1,1)))^k
end

function Et_int(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   return sparse(I,ψb.lng,ψk.lng)
end

function cosα_int(ψb::TPsi,ψk::TPsi,k::Int,q::Int)::SparseMatrixCSC{Float64,Int}
   #if ψb==ψk
      out = spdiagm(k=>fill(0.5,ψb.lng-k),-k=>fill(0.5,ψb.lng-k))
   #else #this won't work correctly. fingers crossed it doesn't come up for now
   #   out = dropzeros!(sparse(δi.(ψb.ms',ψk.ms+k) +δi.(ψb.ms',ψk.ms-k)))
   #end
   l = max(length(ψk.ms[1])*(length(ψk.σ) - 1),1)
   return out = kron(sparse(I,l,l),out)
end
function cosβ_int(ψb::Psi,ψk::Psi,p::Int)::SparseMatrixCSC{Float64,Int}
   if ψb.T==ψk.T
      out = spdiagm(p=>fill(0.5,ψb.T.lng-p),-p=>fill(0.5,ψb.T.lng-p))
   else #this won't work correctly. fingers crossed it doesn't come up for now
      out = dropzeros!(sparse(δi.(ψb.T.ms',ψk.T.ms+p) +δi.(ψb.T.ms',ψk.T.ms-p)))
   end
   l = max(length(ψk.ms[1])*(length(ψk.σ) - 2), 1)
   w = max(length(ψk.ms[1]), 1)
   return out = kron(sparse(I,l,l),out,sparse(I,w,w))
end
function cosγ_int(ψb::Psi,ψk::Psi,p::Int)::SparseMatrixCSC{Float64,Int}
   if ψb.T==ψk.T
      out = spdiagm(p=>fill(0.5,ψb.T.lng-p),-p=>fill(0.5,ψb.T.lng-p))
   else #this won't work correctly. fingers crossed it doesn't come up for now
      out = dropzeros!(sparse(δi.(ψb.T.ms',ψk.T.ms+p) +δi.(ψb.T.ms',ψk.T.ms-p)))
   end
   l = max(length(ψk.ms[1])*(length(ψk.σ) - 3),1)
   w = max(length(ψk.ms[1])*2,1)
   return out = kron(sparse(I,l,l),out,sparse(I,w,w))
end

################################################################################
#####                         Section 3: Dictionary                        #####
################################################################################
function Opsdict()::Dict{String,Function}
out = Dict{String,Function}("E"=>E,"Nz"=>nz, "N2"=>n2,"Np"=>np,"Nm"=>nm,"Npm"=>npm,
    "Nx"=>nx,"iNy"=>iny,"NS"=>ns,"S2"=>s2,"Sz"=>sz,"Sp"=>sp,"Sm"=>sm,"Spm"=>spm,
    "Sx"=>sx,"iSy"=>isy)#,"Pα"=>Pα,"cosα"=>cosα,"Pβ"=>Pβ,"cosβ"=>cosβ,"Pγ"=>Pγ,
    #"cosγ"=>cosγ,"μz"=>μz,"μx"=>μx,"iμy"=>iμy)
end

function Intdict()::Dict{String,Function}
out = Dict{String,Function}("E"=>Et_int,"μz"=>μzf,"μx"=>μxf,"iμy"=>iμyf,"cosα"=>cosα_int,
    "cosβ"=>cosβ_int,"cosγ"=>cosγ_int)
end


#= Graveyard
#function sq_op(j::Real,s::Real,q::Int,ψ::Psi)::SparseMatrixCSC{Float64,Int64}
#   l = ψ.lng
#   out = spzeros(l,l)
#   if s≠zero(s)
#      for a ∈ 1:l, b ∈ 1:l
#         if abs(ψ.N[a] - ψ.N[b])≤q && (q+ψ.K[a]-ψ.K[b])==0
#            @views out[b,a] = skqpart(j,s,q,q,ψ.N[b],ψ.K[b],ψ.N[a],ψ.K[a])
#         end
#      end
#      dropzeros!(out)
#      out .*= nred(s)*powneg1(s+j+q)
#      #The q here is because the rank and component are equal for this function
#   else
#      out[diagind(out)] .+= 1.0
#   end
#   return out
#end
#function sz_op(j::Real,s::Real,ψ::Psi,p::Int)
#   l = ψ.lng
#   out = spzeros(l,l)
#   if s≠zero(s)&&p≠0
#      for a ∈ 1:l, b ∈ a:l
#         if abs(ψ.N[a] - ψ.N[b])≤1 && ψ.K[a]==ψ.K[b]
#            @views out[b,a] = s0part(j,s,ψ.N[b],ψ.K[b],ψ.N[a],ψ.K[a])
#         end
#      end
#      dropzeros!(out)
#      out .*= nred(s)*powneg1(s+j+1)
#      out = Symmetric(out,:L)^p
#   else
#      out[diagind(out)] .= 1.0
#   end
#   return out
#end
#sz(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = sz_op(ψ.J,ψ.S,ψ,p)
################################################################################
#####                         Section 2: Operators                         #####
################################################################################
#E = Op(1.0,tp=[0;0;;])
#Nz = Op(1.0,d=1)
#N2 = Op(1.0,a=1)
#Np = Op(1.0,rp=[1],rf=[np])
#Nm = Op(1.0,rp=[1],rf=[nm])
#Npm = Op(1.0,rp=[1],rf=[npm])
#Nx = Op(1.0,rp=[1],rf=[nx])
#iNy = Op(1.0,rp=[1],rf=[iny])
#
#NS = Op(1.0,b=1)
#S2 = Op(1.0,c=1)
#Sz = Op(1.0,rp=[1],rf=[sz])
#Sp = Op(1.0,rp=[1],rf=[sp])
#Sm = Op(1.0,rp=[1],rf=[sm])
#Spm = Op(1.0,rp=[1],rf=[spm])
#Sx = Op(1.0,rp=[1],rf=[sx])
#iSy = Op(1.0,rp=[1],rf=[isy])
#
#Pα = Op(1.0,tp=[1;0;;])
#cosα = Op(1.0,tp=[0;1;;])
#Pβ = Op(1.0,tp=[0 1; 0 0])
#cosβ = Op(1.0,tp=[0 0; 0 1])
#Pγ = Op(1.0,tp=[0 0 1; 0 0 0])
#cosγ = Op(1.0,tp=[0 0 0; 0 0 1])
#
#μz = Op(1.0,rp=[1],rf=[μzf])
#μx = Op(1.0,rp=[1],rf=[μxf])
#iμy = Op(1.0,rp=[1],rf=[iμyf])

=#