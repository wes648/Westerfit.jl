
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
