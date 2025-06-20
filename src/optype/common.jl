
"""
powneg1(x) takes a number and returns (-1)^x. I realize this is a stupid looking
function to have but it evalutes every so slightly faster than just (-1)^x
"""
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
"""
δi(x,y) is the Kronecker delta that returns an integer type
"""
δi(x::Real,y::Real)::Int = x==y

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
function srprep2(J::Number,S::Number)
   ns = Δlist2(J,S)
   #nd = 2 .* Int.(ns) .+ 1
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


function adiagin(a::AbstractMatrix)::StepRange{Int,Int}
   @assert size(a,1)==size(a,2)
   l = size(a,1)
   return l:(l-1):(l*(l-1)+1)
end

"""
This builds the rotational Wang Transformation matrix for every n in Δlist(j,s).
"""
function ur(j::Float64,s::Float64)::SparseMatrixCSC{Float64, Int64}
   ns, ni, jsd = srprep2(j,s)
   out = spzeros(jsd,jsd)
   for i in 1:length(ns)
      out[ni[i], ni[i]] = ur(ns[i])
   end
   return out
end
function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   out = Diagonal(append!(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(append!(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end
#function ur(n::Int)::SparseMatrixCSC{Float64, Int}
#   return sparse(I,2n+1,2n+1)
#end
function ur2(n::Int)::SparseMatrixCSC{Float64, Int}
#this version is radically faster above N = 20 but much slower
   out = spdiagm(append!(fill(-√.5,n), 0.0, fill(√.5,n)))
   out[adiagin(out)] .= append!(fill(√.5,n), 1.0, fill(√.5,n))
   return sparse(out)
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

Δlist2(J::Real,S::Real)::UnitRange{Int} = Int(abs(J-S)):Int(J+S)
kgen(ns::UnitRange{Int})::Vector{UnitRange{Int}} = [-n:n for n ∈ ns]

#tplus!(a::Diagonal)::SparseMatrixCSC{Float64, Int} = sparse(a)
#tplus!(a::Array{Float64,2})::Array{Float64,2} = hermitianpart!(2a)
function tplus!(a::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   a .+= permutedims(a)
end
function tplus!(a::SparseMatrixCSC{ComplexF64,Int})::SparseMatrixCSC{ComplexF64,Int}
   a .+= permutedims(conj(a))
end
