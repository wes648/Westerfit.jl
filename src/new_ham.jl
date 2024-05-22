

eh(x::Number)::Float64 = √(x*(x+1))
□rt(x::Number)::Float64 =√(x*(x>zero(x)))
fh(x::Number,y::Number)::Float64 = □rt((x-y)*(x+y+1))
jnred(j::Number,n::Number)::Float64 = √((2*j+1)*(2*n+1))
nred(n::Number)::Float64 = √(n*(n+1)*(2*n+1))
powneg1(k::Number)::Int = isodd(k) ? -1 : 1
Δlist(J,S)::Array{Int64} = collect(Int(abs(J-S)):Int(J+S))

function srprep(J,S)
   ns = Δlist(J,S)
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
function qngen(j,s)
   ns, nd, ni, jsd = srprep(j,s)
   out = zeros(Int,jsd,2)
   for i in 1:length(ns)
      out[ni[i,1]:ni[i,2],1] .= ns[i]
      out[ni[i,1]:ni[i,2],2] = collect(Int,-ns[i]:ns[i])
   end
   #[n k]
   return out
end

function nnss_check(a,b)::Int
   a = a*iseven(a) + (a-1)*isodd(a)
   b = b*iseven(b) + (b-1)*isodd(b)
   return min(a,b)
end
ns_el(j,s,p,n) = (0.5*eh(j) - eh(n) - eh(s))^p
function nnss_op(j,s,a,b,qns)
   c = nnss_check(a,b)
   a -= c
   b -= c
   @views out = eh.(qns[:,1]).^a .* ns_el.(j,s,c,qns[:,1]) .* eh(s)^b
   return Diagonal(out)
end

nz_op(qns,p) = @views out = Diagonal(qns[:,2].^p)

function np_op(j,s,qns,p)
   ns = qns[1+p:end,1]
   ks = qns[1+p:end,2]
   out = fh.(ns,ks.-1)
   for o in 2:p
      out .*= fh.(ns,qns[1+p:end,2].-o)
   end
   out = spdiagm(-p=>out)
   return out
end
npm_op(j,s,qns,p) = Symmetric(np_op(j,s,qns,p),:L)

function sqpart(j,s,q,bqn,kqn)::Float64
   nb = bqn[1]
   kb = bqn[2]
   nk = kqn[1]
   kk = kqn[2]
   return wig3j(nb,1,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,1)*jnred(nb,nk)*powneg1(-kb)
end
function sq_op(j,s,q,qns)#::SparseMatrixCSC{Float64, Int64}
   l = size(qns,1)
   out = spzeros(l,l)
   kcol = spzeros(Int,l,l)
   if s != zero(s)
      for a ∈ 1:l, b ∈ 1:l
#         if abs(qns[a,1]-qns[b,1])≤1 && (q+qns[a,2]-qns[b,2])==0
            println("a = $a, b =$b")
            println("N' = $(qns[b,1]), N = $(qns[a,1]), K' = $(qns[b,2]), K = $(qns[a,2])")
            @views out[b,a] = sqpart(j,s,q,qns[b,:],qns[a,:])
            kcol[b,a] = qns[a,2]
#         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1)
   else
      out[diagind(out)] .+= 1.0
   end
   return out, kcol
end


      #if q==0
      #   out = Symmetric(out,:L)
      #else
      #   out += transpose(out)
      #end
