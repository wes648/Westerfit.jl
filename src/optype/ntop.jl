
"""
Returns the number of σ states for the value of nfold[1]
"""
function σcount(nfold::Real)::Int
   if isodd(nfold)
      out = ceil(Int,0.5*nfold)
   elseif iseven(nfold)&&nfold≠0
      out = floor(Int,nfold/4) + 1
   else #nfold == 0
      out = 1
   end
   return out
end
"""
The counting of number of sigma states is different for additional rotors
"""
function σcount_post(nfold::Real)::Int
   if nfold ≤ 2
      out = 1
   elseif nfold==3
      out = 3
   else
      @warn "secondary tops must be nfold ≤ 3 for now"
   end
   return out
end
function σcount(nfold::Vector{Int})::Int
   out = σcount(nfold[1])
   for i ∈ 2:length(nfold)
      out += σcount_post(nfold[i])
   end
   return out
end

"""
Generates all the sigmas for the 2 rotor problem based on both sym folds
"""
function σgen(nf1::Int,nf2::Int)::Array{Int,2}
   σcnt1 = σcount(nf1)
   σcnt2 = σcount(nf2)-1
   out = zeros(Int,2,0)
   for i in 1:σcnt1
      σ1 = i-1
      if σ1==0
         σ2 = collect(0:σcnt2)'
      else
         σ2 = collect(σcnt2:-1:-σcnt2)'
      end#if
      out = hcat(out, [fill(σ1,length(σ2))'; σ2])
   end#for
   return out
end#function
"""
Generates all the σs for a vector of rotor folds
"""
function σgen(nf::Array{Int})::Array{Int,2}
   if length(nf) > 1
   old = σgen(nf[1],nf[2])
   for i in 3:length(nf)
      σlsti = nf[i] - iseven(nf[i])
      new = vcat(old,fill(0,size(old,2))')
      for j in 2:σlsti
         σ = powneg1(j)*floor(Int,j/2)
         part = vcat(old[:,1+(j==σlsti):end],fill(σ,size(old[:,1+(j==σlsti):end],2))')
         new = hcat(new,part)
      end#for j
      old = new
   end#for i
   else
      old = σgen(nf[1])
   end
   return old
end
"""
Generates the σs for 1 rotor
"""
function σgen(nf::Int)::Array{Int}
   return collect(0:σcount(nf[1])-1)'
end
"""
generates list of σs for a later rotor
"""
function σgen_post(nf::Int)::Array{Int}
   if nf ≤ 2
      out = [0]
   elseif nf==3
      out = [1;0;-1]
   else
      @warn "secondary tops must be nfold ≤ 3 for now"
   end
   return out
end
"""
Generates the σs for the nth top
"""
function nth_σgen(nfold::Int,tid::Int)::Array{Int}
   if isone(tid)
      return σgen(nfold)
   else
      return σgen_post(nfold)
   end
end
"""
Generates the number of σs for the nth top
"""
function nth_σcount(nfold::Int,tid::Int)::Int
   if isone(tid)
      return σcount(nfold)
   else
      return σcount_post(nfold)
   end
end

"""
Returns the list of ms for a gien nfold, basis size, and σ
"""
function msgen(nfold::Int,mcalc::Int,σ::Int)::StepRange{Int,Int}
   if iszero(nfold)
      marray = 0:0
   elseif isodd(nfold)
      lim = mcalc*nfold
      marray = (-lim+σ):nfold:(lim+σ)
   else #iseven(nfold)
      lim = mcalc*nfold
      lim = floor(Int,lim/2)
      marray = -lim:floor(Int,nfold/2):lim
      marray = (-lim+σ):nfold:(lim+σ)
   end
   return marray
end
#msgen(nfold::Vector{Int},mcalc::Int,σ::Int)::StepRange{Int,Int} = msgen(nfold[1],mcalc,σ)

function msgen(nf::Array{Int},mcalc::Int,σs::Array{Int})::Vector{StepRange{Int,Int}}
   out = Vector{StepRange{Int,Int}}(undef,length(nf))
   for j in eachindex(nf)
      out[j] = msgen(nf[j],mcalc,σs[j])
   end
   return out
end

"""
σ2ind returns the index of the σ states for the ith top of n-fold symmetry.
Args: σ::Int, tid::Int, nfold::Int
"""
function σ2ind(σ::Int,tid::Int,nfold::Int)::Int
   if isone(tid)
      return σ+1
   elseif tid > 1
      return σcount(nfold) - σ
   else
      @warn "why is tid zero?"
   end
end
function σ2ind(ψ::TTPsi,tid::Int)::Int
   if isone(tid)
      return ψ.σs[1]+1
   elseif tid > 1
      return σcount(ψ.nfs[tid]) - ψ.σs[tid]
   else
      @warn "why is tid zero?"
   end
end

"""
I don't understand what or why this function exists
it appears to return 
"""
function nσfinder(tid::Int,σ::Int,nfold::Vector{Int})::Int
   if isone(tid)
      out = σ+1
   elseif !iszero(tid)
      out = σcount(nfold[1])
      out += sum(nfold[2:tid-1])
      out += σ+2
   else
      @warn "top id is zero. this will crash"
   end
   return out
end
