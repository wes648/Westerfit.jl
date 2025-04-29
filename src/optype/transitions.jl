

struct μFunc
   μv::Float64
   fr::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,Int}}
   ft::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{TPsi,Int}}
   k::Int
   q::Int
   function μFunc(μv::Float64,fr::Function,ft::Function,k::Int,q::Int)
      new(μv,fr,ft,k,q)
   end
end

function eval_μop_r(op::μFunc,ψb::RPsi,ψk::RPsi)::SparseMatrixCSC{Float64,Int} 
   op.fr(ψb,ψk,op.k,op.q)
end
function eval_μop_t(op::μFunc,ψb::TPsi,ψk::TPsi)::SparseMatrixCSC{Float64,Int} 
   op.ft(ψb,ψk,op.k,op.q)
end

function int_enact(μp,μf,ψb,ψk,bvec,kvec)::SparseMatrixCSC{Float64,Int}
   out = kron(μp,eval_μop_t(μf,ψb.T,ψk.T),eval_μop_r(μf,ψb.R,ψk.R))
end

function bld_μmat_1s(inthres,μp,μf,ψb,ψk,bvec,kvec)::SparseMatrixCSC{Float64,Int}
   out = int_enact(μp[1],μf[1],ψb,ψk,bvec,kvec)
   for i ∈ 2:length(μf)
      out += int_enact(μp[i],μf[i],ψb,ψk,bvec,kvec)
   end
   droptol!(out,inthres)
   return out
end

function jbjklister(jmin,jmax,mΔj)
   out = zeros(0,2)
   for j in jmin:(jmax-1)
      out = vcat(out,[collect(j:(mΔj+j)) fill(j,mΔj+1)])
   end
   out = vcat(out,[collect(jmax:(mΔj+jmax-1)) fill(jmax,mΔj)])
   return out
end

function intcalc()
   s = ctrl["S"]
   vtm = ctrl["vtmax"]
   nf = ctrl["NFOLD"]
   mΔj = 1
   jbjk = jbjklister(0.5*iseven(Int(2*s+1)),jmax,mΔj)

   for i in 1:size(jbjk,1)
      jb = jbjk[i,1]
      jk = jbjk[i,2]
      #println("jb=$jb, jk=$jk")
      #needs to be corrected for the vt scheme
      kinds = jvlinds(jk,s,vtm)
      binds = jvlinds(jb,s,vtm)
      #filter vecs
      kvecs = vecs[1:Int(2*jk+1)*cmsd,kinds]
      bvecs = vecs[1:Int(2*jb+1)*rmsd,binds]
      ψb = Psi( RPsi(jb,s), TPsi(nf,σ,mc) )
      ψk = Psi( RPsi(jb,s), TPsi(nf,σ,mc) )
      #calculate intensities
      μs = intmat(ctrl["INTthres"],μp,μf,ψb,ψk,bvecs,kvecs)

      if (jb==jk)#&&(σr==σc)
         μs = sparse(UpperTriangular(μs))
      end
      rind, cind, tints = findnz(μs)
      for l in 1:length(tints)
         b = binds[rind[l]]
         k = kinds[cind[l]]
         ints[b,k] = tints[l]
      end
      #ints[binds[cind],kinds[rind]] = tints
   end
   return ints
end

function tracalc()
   ints = intcalc()
   rinds, cinds, = findnz(ints)
   for i in 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = rvals[r] - cvals[c]
      frqs[r,c] = ν*(ctrl["νmin"]≤abs(ν*0.001)≤ctrl["νmax"])
   end
   droptol!(frqs,1e-3) #drops all frequencies below 1kHz
   #assign Elow & qns plus sign correction
   rinds, cinds, νs = findnz(frqs)
   len = nnz(frqs) 
   outfrq = zeros(len,4)
   outqns = zeros(Int,len,12)
   for i in 1:len
      ν = νs[i]
      r = rinds[i]
      c = cinds[i]
      if ν > 0.0
         outfrq[i,1] = ν
         outfrq[i,3] = cvals[c] / csl
         outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-cvals[c]/kbT) - exp(-rvals[r]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σc+1)*100
         outqns[i,1:6] = rqns[r,:]
         outqns[i,7:12] = cqns[c,:]
      elseif ν < 0.0
         outfrq[i,1] = -ν
         outfrq[i,3] = rvals[r] /csl
         outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-rvals[r]/kbT) - exp(-cvals[c]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*(2*s+1)*(σr+1)*100
         outqns[i,1:6] = cqns[c,:]
         outqns[i,7:12] = rqns[r,:]
      end
   end
   #println(outfrq)
   #println(outqns)
   return outfrq, outqns
end

end

################################################################################
###     Parametric Typing Testing
#=
using FunctionWrappers
import FunctionWrappers: FunctionWrapper

struct StabOp{T<:Number}
   f::FunctionWrapper{T, Tuple{Float64,Int}}
   p::Int
   function StabOp(T::Type,f::Function,p::Int)
      new{T}(f,p)
   end
end

function eval_op(op::StabOp{T},arg)::T where {T<:Number}
   op.f(arg, op.p)
end

e(x::Float64,p::Int) = x^p
f(x::Float64,p::Int) = im * (x^p)
g(x::Float64,p::Int) = cos(p*x)
h(x::Float64,p::Int) = im*sin(p*x)

E = StabOp(Float64,e,2)
F = StabOp(ComplexF64,f,3)
G = StabOp(Float64,g,4)
H = StabOp(ComplexF64,h,5)
=#