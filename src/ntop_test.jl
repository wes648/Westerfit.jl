using LinearAlgebra,SparseArrays
function pa_ntop(nf,ms,tid,p)::Diagonal{Float64}
   out = ms[:,tid] .^ p
   out = ntopset!(nf,ms,tid,out)
   out = Diagonal(out)
   return out
end
function cos_ntop(nf,ms::Array{Int},tid,p::Int)::SparseMatrixCSC{Float64, Int64}
   if p==0
      out = I(size(ms,1))
   else
      out = fill(0.5, size(ms,1)-p)
      out = spdiagm(-p=>out, p=>out)
      out = ntopset!(nf,ms,tid,out)
   end
   return out
end

function htest_1stg(mcalc)
   Fa = 5.1
   Fb = 5.2
   Fab = 0.
   V3a = 100.
   V3b = 300.
   Vab = 0.
   nf = [3;3]
   σs = σgen(3,3)
   for i in 1:size(σs,2)
      ms = msgen_indef(nf,mcalc,σs[:,i])
      H = Fa*pa_ntop(nf,ms,1,2) + Fb*pa_ntop(nf,ms,2,2) + Fab*pa_ntop(nf,ms,1,1)*pa_ntop(nf,ms,2,1)
      H += -0.5*(V3a*cos_ntop(nf,ms,1,1) + V3b*cos_ntop(nf,ms,2,1)) + Vab*cos_ntop(nf,ms,1,1)*cos_ntop(nf,ms,2,1)
      H += 0.5*(V3a+V3b)*I
      @show σs[:,i]
      H = Symmetric(Matrix(H))
      @show eigvals(H)
   end
end