"""
Welcome to be base Operators used in westerfit! 
This file is divided into 3 sections and to add a new Operator, all 3 must be
  edited.
Section 1: Functions. These act on the wavefuntions to produce the matrices and
must use the following type annotations:
 op(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64}
This should ensure consistency with the program.
Stage 0 operators (those used in the intensity calculations) have the form of:
 op(ψb::Psi,ψk::Psi,p::Int)::SparseMatrixCSC{Float64, Int64}
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
ns_el(j,s,p,n)::Float64 = (0.5*eh(j) - eh(n) - eh(s))^p

#These functions are cooked into enact_init as they are purely diagonal
#they are N^2a (NS)^b S^c Nz^d in this order
function n2(ns::UnitRange{Int},p::Int)::Vector{Float64}
   reduce(vcat, [fill(eh(n)^p,2n+1) for n ∈ ns])
end
function n2(v::Float64,ns::UnitRange{Int},p::Int)::Vector{Float64}
   reduce(vcat, [fill(x*eh(n)^p,2n+1) for n ∈ ns])
end
function ns_el2(j,s,ns::UnitRange,p::Int)::Vector{Float64}
   reduce(vcat, fill[(eh(j) - eh(n) - eh(s))^p, 2n+1 for n ∈ ns])
end
nz(ks::UnitRange{Int},p::Int)::Vector{Float64} = reduce(vcat, ks).^p

function np(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64}
   ns = ψ.N[1+p:end]
   part = ones(length(ns))
   if p ≤ ψ.lng && p ≠ 0
      ks = ψ.K[1+p:end]
      part = ones(length(ks))
      for o in 1:p
         part .*= fh.(ns,ks.-o)
      end
   #end#original if
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

function npv2(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64}
   ns = ψ.N[1+p:end]
   part = ones(length(ns))
   if p ≤ ψ.lng && p ≠ 0
      ks = ψ.K[1+p:end]
      part = ones(length(ks))
      for n in ns #this won't work correctly but might be closer to the ideal
         part = prod(fh.(n,collect(ks) .- collect(1:o)),dims=2)
      end
      out = spdiagm(-p=>part)
   elseif p > ψ.lng && p ≠ 0
      out = spzeros(ψ.lng,ψ.lng)
   else
      out = spdiagm(ones(ψ.lng))
   end
   return out
end


nm(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64} = permutedims(np(ψ,p))
npm(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64} = tplus!(np(ψ,p))
nx(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64} = tplus!(0.5*np(ψ,p))

function iny(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int}
   if p≠0
      out = np_op(ψ,1-δi(p,0))
      out .-= permutedims(out)
      return dropzeros!(out)
   else
      return spdiagm(ones(ψ.lng))
   end
end
#iny(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64} = 

function skqpart(j,s,k,q,nb,kb,nk,kk)::Float64
   return wig3j(nb,k,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,k)*jnred(nb,nk)*powneg1(-kb)
end
function s0part(j,s,nb,kb,nk,kk)::Float64
   return skqpart(j,s,1,0,nb,kb,nk,kk)::Float64
end
function sz_op(j::Real,s::Real,ψ::Psi,p::Int)
   l = ψ.lng
   out = spzeros(l,l)
   if s≠zero(s)&&p≠0
      for a ∈ 1:l, b ∈ a:l
         if abs(ψ.N[a] - ψ.N[b])≤1 && ψ.K[a]==ψ.K[b]
            @views out[b,a] = s0part(j,s,ψ.N[b],ψ.K[b],ψ.N[a],ψ.K[a])
         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1)
      out = Symmetric(out,:L)^p
   else
      out[diagind(out)] .= 1.0
   end
   return out
end
sz(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = sz_op(ψ.J,ψ.S,ψ,p)

function sq_op(j::Real,s::Real,q::Int,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   l = ψ.lng
   out = spzeros(l,l)
   if s≠zero(s)
      for a ∈ 1:l, b ∈ 1:l
         if abs(ψ.N[a] - ψ.N[b])≤q && (q+ψ.K[a]-ψ.K[b])==0
            @views out[b,a] = skqpart(j,s,q,q,ψ.N[b],ψ.K[b],ψ.N[a],ψ.K[a])
         end
      end
      drOpzeros!(out)
      out .*= nred(s)*powneg1(s+j+q)
      #The q here is because the rank and component are equal for this function
   else
      out[diagind(out)] .+= 1.0
   end
   return out
end
sp(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = p*sq_op(ψ.J,ψ.S, p,ψ)
sm(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = p*sq_op(ψ.J,ψ.S,-p,ψ)
spm(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int} = tplus!(p*sq_op(ψ.J,ψ.S, p,ψ))
function sx(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int}
   out  = sq_op(ψ.J,ψ.S, 1,ψ)
   out -= sq_op(ψ.J,ψ.S,-1,ψ)
   out *= -√.5
   out ^= p
   return out
end
function isy(ψ::Psi,p::Int)::SparseMatrixCSC{Float64,Int}
   out  = sq_op(ψ.J,ψ.S, 1,ψ)
   out += sq_op(ψ.J,ψ.S,-1,ψ)
   out *= -√.5
   out ^= p
   return out
end

function μzf(ψb::Psi,ψk::Psi,p)::SparseMatrixCSC{Float64,Int}
   spdiagm(ones(ψk.lng))
end
function μxf(ψb::Psi,ψk::Psi,p)::SparseMatrixCSC{Float64,Int}
   spdiagm(ones(ψk.lng))
end
function iμyf(ψb::Psi,ψk::Psi,p)::SparseMatrixCSC{Float64,Int}
   spdiagm(ones(ψk.lng))
end


################################################################################
#####                         Section 2: Operators                         #####
################################################################################
E = Op(1.0)
Nz = Op(1.0,d=1)
N2 = Op(1.0,a=1)
Np = Op(1.0,rp=[1],rf=[np])
Nm = Op(1.0,rp=[1],rf=[nm])
Npm = Op(1.0,rp=[1],rf=[npm])
Nx = Op(1.0,rp=[1],rf=[nx])
iNy = Op(1.0,rp=[1],rf=[iny])

NS = Op(1.0,b=1)
S2 = Op(1.0,c=1)
Sz = Op(1.0,rp=[1],rf=[sz])
Sp = Op(1.0,rp=[1],rf=[sp])
Sm = Op(1.0,rp=[1],rf=[sm])
Spm = Op(1.0,rp=[1],rf=[spm])
Sx = Op(1.0,rp=[1],rf=[sx])
iSy = Op(1.0,rp=[1],rf=[isy])

Pα = Op(1.0,tp=[1;0;;])
cosα = Op(1.0,tp=[0;1;;])
Pβ = Op(1.0,tp=[0 1; 0 0])
cosβ = Op(1.0,tp=[0 0; 0 1])
Pγ = Op(1.0,tp=[0 0 1; 0 0 0])
cosγ = Op(1.0,tp=[0 0 0; 0 0 1])

μz = Op(1.0,rp=[1],rf=[μzf])
μx = Op(1.0,rp=[1],rf=[μxf])
iμy = Op(1.0,rp=[1],rf=[iμyf])

################################################################################
#####                         Section 3: Dictionary                        #####
################################################################################
function Opsdict()::Dict{String,Op}
out = Dict{String,Op}("E"=>E,"Nz"=>Nz, "N2"=>N2,"Np"=>Np,"Nm"=>Np,"Npm"=>Npm,
    "Nx"=>Nx,"iNy"=>iNy,"NS"=>NS,"S2"=>S2,"Sz"=>Sz,"Sp"=>Sp,"Sm"=>Sm,"Spm"=>Spm,
    "Sx"=>Sx,"iSy"=>iSy,"Pα"=>Pα,"cosα"=>cosα,"Pβ"=>Pβ,"cosβ"=>cosβ,"Pγ"=>Pγ,
    "cosγ"=>cosγ,"μz"=>μz,"μx"=>μx,"iμy"=>iμy)
end