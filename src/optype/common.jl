
"""
powneg1(x) takes a number and returns (-1)^x. I realize this is a stupid looking
function to have but it evalutes every so slightly faster than just (-1)^x
"""
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
"""
δi(x,y) is the Kronecker delta that returns an integer type
"""
δi(x::Real,y::Real)::Int = x==y

Δlist2(J::Real,S::Real)::UnitRange{Int} = Int(abs(J-S)):Int(J+S)

dgen(x::Real)::Int = Int(2*x) + 1

sand(a::AbstractArray,b::AbstractArray) = b' * a * b

function srprep(J,S)
   ns = Δlist2(J,S)
   nd = 2 .* Int.(ns) .+ 1
   ni = ones(Int, length(ns),2)
   ni[1,2] = nd[1]
   for i in 2:length(ns)
      ni[i,1] = ni[i-1,2] + 1
      ni[i,2] = ni[i,1] + nd[i] - 1
   end
   jd = Int((2.0*S+1.0)*(2.0*J+1.0))
   return ns, nd, ni, jd
end
function srprep2(J::Real,S::Real)
   ns = Δlist2(J,S)
   ni = nindsgen(ns)
   jd = Int((2.0*S+1.0)*(2.0*J+1.0))
   return ns, ni, jd
end

function nindsgen(ns::UnitRange{Int})::Vector{UnitRange{Int}}
   ni = Vector{UnitRange{Int}}(undef,length(ns))
   ni[1] = 1:2*ns[1]+1
   for i in 2:length(ns)
      ni[i] = (ni[i-1][end]+1):(ni[i-1][end]+2ns[i]+1)
   end
   return ni
end

"""
Returns a steprange of the indices of along the anti-diagonal of the input matrix
"""
function adiagin(a::AbstractMatrix)::StepRange{Int,Int}
   @assert size(a,1)==size(a,2)
   l = size(a,1)
   return l:(l-1):(l*(l-1)+1)
end

"""
This builds the rotational Wang Transformation matrix for every n in Δlist(j,s).
"""
function ur(j::Real,s::Real)::SparseMatrixCSC{Float64, Int64}
   if !iszero(s)
      ns, ni, jsd = srprep2(j,s)
      out = spzeros(jsd,jsd)
      for i in 1:length(ns)
         out[ni[i], ni[i]] = ur(ns[i])
      end
      return out
   else
      return ur(Int(j))
   end
end
function ur_old(n::Int)::SparseMatrixCSC{Float64, Int}
   out = Diagonal(append!(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(append!(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end
#function ur(n::Int)::SparseMatrixCSC{Float64, Int}
#   return sparse(I,2n+1,2n+1)
#end
function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   if !iszero(n)
   out = spzeros(1:2n+1,1:2n+1)
   out[diagind(out)[1:n]] .= -√.5
   out[diagind(out)[n+2:end]] .= √.5
   out[adiagin(out)] .= √.5
   out[n+1,n+1] = 1.0
   return sparse(out)
   else
      return sparse(1.0I,1,1)
   end
end

#ur(n::Int)::SparseMatrixCSC{Float64, Int} = n ≤ 20 ? ur1(n) : ur2(n)
#ur(n::Int)::SparseMatrixCSC{Float64, Int} = ur(n)

function jvdest2(j::Float64,s::Float64,vtm::Int)::UnitRange{Int}
"""
This returns a unit range spanning from the first to the final indices for a 
   certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (vtm+1)*(2*s+1)*2*sum(collect((0.5*isodd(2*s)):(j-1)) .+0.5))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*2*sum(collect((0.5*isodd(2*s)):j) .+0.5))
   return snd:fnd
end

function qnlab(j,s,vtm,σ)::Array{Int,2}
   nlist = Δlist2(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   vd = Int(vtm+1)
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out[:,2] = abs.(out[:,2])
   out = kron(ones(Int,vd),out)
   vtrray = kron(collect(0:vtm) ,ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
   return out
end
"""
Determines the value of Kc based on the value of N and |Kₐ|
"""
function k2kc(n,k)
   ka = abs(k)
   if k < 0
      kc = n - ka + 1 - isodd(n + k)
   elseif k == zero(k)
      kc = n
   else
      kc = n - ka + isodd(n + k)
   end
   return kc
end

function bigqngen(l,jlist,s,vtm,σs)
   σcnt = maximum(size(σs))
   out = zeros(Int,l,6,σcnt)
   for j ∈ jlist[:,1]
      dest = jvdest2(0.5*j,s,vtm)
      for σ in 1:σcnt
      out[dest,:,σ] = qnlab(0.5*j,s,vtm,σ-1)
   end;end
   return out
end

kgen(ns::UnitRange{Int})::Vector{UnitRange{Int}} = [-n:n for n ∈ ns]

#tplus!(a::Diagonal)::SparseMatrixCSC{Float64, Int} = sparse(a)
#tplus!(a::Array{Float64,2})::Array{Float64,2} = hermitianpart!(2a)
function tplus!(a::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   a .+= permutedims(a)
end
function tplus!(a::SparseMatrixCSC{ComplexF64,Int})::SparseMatrixCSC{ComplexF64,Int}
   a .+= permutedims(conj(a))
end


#indexes for ntop operators
#ti is the top index, nt is the number of tops
ffind(ti::Int,nt::Int)::Int = 11 + ti
rzind(ti::Int,nt::Int)::Int = 11 + ti +   nt
rxind(ti::Int,nt::Int)::Int = 11 + ti + 2*nt
vnind(ti::Int,nt::Int)::Int = 11 + ti + 3*nt
ezind(ti::Int,nt::Int)::Int = 11 + ti + 4*nt
exind(ti::Int,nt::Int)::Int = 11 + ti + 5*nt

"""
Applies Kronecker products with identity matrices in order to properly resize the ith one top matrix.
"""
function torsetter!(ψ::TPsi,i::Int,out)
   lnf = length(ψ.nf)
   lbk = ψ.lb
   if lnf > 1
      out = kron( sparse(I, lbk*(lnf-i), lbk*(lnf-i)), 
                  out, 
                  sparse(I, lbk*(i-1), lbk*(i-1)) )
   end
   return out
end

"""
A simple wrapper for making the matrix Symmetric if Real or Hermitian if complex then diagonalize.
Makes life easier and helps me keep the fully real setup for C_s while permitting the flexibility for C_1
"""
function diagwrap(H::AbstractArray)::Eigen
   if eltype(H)<:Real #isreal(eltype(H))
      return eigen!(Symmetric(Matrix(H), :L))
   else
      return eigen!(Hermitian(Matrix(H), :L))
   end
end
