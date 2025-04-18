
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

function σgen_indef(nf::Array{Int})::Array{Int,2}
   old = σgen(nf[1],nf[2])
   for i in 3:length(nf)
      σlsti = nf[i] - iseven(nf[i])
      new = vcat(old,fill(0,size(old,2))')
      for j in 2:σlsti
         σ = Int(cospi(j)*floor(j/2))
         part = vcat(old[:,1+(j==σlsti):end],fill(σ,size(old[:,1+(j==σlsti):end],2))')
         new = hcat(new,part)
      end#for j
      old = new
   end#for i
   return old
end
function σgen_indef(nf::Int)::Array{Int}
   return collect(0:σcount(nf[1])-1)'
end

function msgen(nfold::Int,mcalc::Int,σ::Int)::StepRange{Int,Int}
   if iszero(nfold)
      marray = 0:0
   elseif isodd(nfold)
      lim = mcalc*nfold
      marray = (-lim+σ):nfold:(lim+σ)
   else #iseven(nfold)
      lim = floor(Int,lim/2)
      marray = -lim:floor(Int,nfold/2):lim
      marray = (-lim+σ):nfold:(lim+σ)
   end
   return marray
end

function msgen_indef(nf::Array{Int},mcalc::Int,σs::Array{Int})::Vector{StepRange{Int,Int}}
   out = Vector{StepRange{Int,Int}}(undef,length(nf))
   for j in eachindex(nf)
      out[j] = msgen(nf[j],mcalc,σs[j])
   end
   return out
end
