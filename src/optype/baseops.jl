"""
Welcome to be base Operators used in westerfit! 
This file is divided into 3 sections and to add a new Operator, all 3 must be
  edited.
Section 1: Functions. These act on the wavefuntions to produce the matrices and
must only take the arguments of a Psi and a power

Section 2: Operators. These are the actual objects used by the code. Look at the
type file to see the structure but the general set up should be:
   Nam = Op(1.0,rp=[1],rf=[function])
where Nam is the name that will be read from the input file and function is the
new Operator function. We are defining it as being raised to the 1st power with
a default parameter value of 1.0

Section 3: Dictionary. This is the list of Operators read by the code in order
to interpret the input file. It is wrapped inside a function so that it doesn't 
stay in memory after being used in the inpu reading.
Each entry should just be "Nam"=>Nam
"""

################################################################################
#####                         Section 1: Functions                         #####
################################################################################
eh(x::Real)::Float64 = x*(x+1)
□rt(x::Real)::Float64 =√(x*(x>zero(x)))
fh(x::Real,y::Real)::Float64 = □rt((x-y)*(x+y+1))
ns_el(j,s,p,n)::Float64 = (0.5*eh(j) - eh(n) - eh(s))^p

#nz(ψ::Psi,p::Int)::Diagonal{Float64} = Diagonal(ψ.K .^p)
#nt(ψ::Psi,p::Int)::Diagonal{Float64} = Diagonal(eh.(ψ.N) .^p)
#ns(ψ::Psi,p::Int)::Diagonal{Float64} = Diagonal(eh(ψ.J) - eh(ψ.S) - eh.(ψ.N))^p
#st(ψ::Psi,p::Int)::Diagonal{Float64} = Diagonal(eh.(ψ.N) .^p)
function np(ψ,p::Int)::SparseMatrixCSC{Float64, Int64}
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
nm(ψ::Psi,p::Int)::SparseMatrixCSC{Float64, Int64} = permutedims(np(ψ,p))
npm(ψ::Psi,p::Int) = tplus!(np(ψ,p))

function iny_op(qns::Matrix{Int64},p::Int)::SparseMatrixCSC{Float64,Int}
   if p≠0
      out = np_op(qns,1-δi(p,0))
      out .-= permutedims(out)
      return dropzeros!(out)
   else
      return spdiagm(ones(ψ.lng))
   end
end
function sqpart(j,s,q,nb,kb,nk,kk)::Float64
   return wig3j(nb,1,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,1)*jnred(nb,nk)*powneg1(-kb)
end
function sz_op(j::Real,s::Real,ψ::Psi,p::Int)
   l = ψ.lng
   out = spzeros(l,l)
   if s≠zero(s)&&p≠0
      for a ∈ 1:l, b ∈ a:l
         if abs(ψ.N[a] - ψ.N[b])≤1 && ψ.K[a]==ψ.K[b]
            @views out[b,a] = sqpart(j,s,0,ψ.N[b],ψ.K[b],ψ.N[a],ψ.K[a])
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
sz(ψ,p)::SparseMatrixCSC{Float64,Int} = sz_op(ψ.J,ψ.S,ψ,p)

function sq_op(j,s,q,qns)::SparseMatrixCSC{Float64, Int64}
   l = size(qns,1)
   out = spzeros(l,l)
   if s≠zero(s)
      for a ∈ 1:l, b ∈ 1:l
         if abs(qns[a,1]-qns[b,1])≤1 && (q+qns[a,2]-qns[b,2])==0
            @views out[b,a] = sqpart(j,s,q,qns[b,:],qns[a,:])
         end
      end
      drOpzeros!(out)
      out .*= nred(s)*powneg1(s+j+1+δ(1,q))*√2
      #the √2 is included to convert from spherical to cylinderical 
   else
      out[diagind(out)] .+= 1.0
   end
   return out
end

function sp_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   #if p≠0
      return sq_op(j,s,1,qns)^p
   #else
   #   return spdiagm(ones(size(qns,1)))
   #end
end
function sm_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   return sq_op(j,s,-1,qns)^p
end
function spm_op(j::Real,s::Real,qns::Array{Int,2},p::Int
         )::SparseMatrixCSC{Float64, Int64} 
   if p≠0
      return tplus!(sp_op(j,s,qns,p))# + sm_op(j,s,qns,p)
   else
      return spdiagm(ones(size(qns,1)))
   end
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

μz = Op(1.0,rp=[1],rf=[μz])
μx = Op(1.0,rp=[1],rf=[μx])
iμy = Op(1.0,rp=[1],rf=[iμy])

################################################################################
#####                         Section 3: Dictionary                        #####
################################################################################
function Opsdict()::Dict{String,Op}
out = Dict{String,Op}("E"=>E,"Nz"=>Nz, "N2"=>N2,"Np"=>Np,"Nm"=>Np,"Npm"=>Npm,
    "Nx"=>Nx,"iNy"=>iNy,"NS"=>NS,"S2"=>S2,"Sz"=>Sz,"Sp"=>Sp,"Sm"=>Sm,"Spm"=>Spm,
    "Sx"=>Sx,"iSy"=>iSy,"Pα"=>Pα,"cosα"=>cosα,"Pβ"=>Pβ,"cosβ"=>cosβ,"Pγ"=>Pγ,
    "cosγ"=>cosγ,"μz"=>μz,"μx"=>μx,"iμy"=>iμy)
end