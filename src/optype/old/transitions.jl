

struct μFuncR
   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{RPsi,RPsi,Int,Int}}
   k::Int
   q::Int
   function μFuncR(f::Function,k::Int,q::Int)
      new(f,k,q)
   end
end
struct μFuncT
   f::FunctionWrapper{SparseMatrixCSC{Float64,Int}, Tuple{TPsi,TPsi,Int,Int}}
   k::Int
   q::Int
   function μFuncT(f::Function,k::Int,q::Int)
      new(f,k,q)
   end
end

struct μOb
   v::Float64
   fr::μFuncR
   ft::μFuncT
   function μOb(μv::Float64,fr::Function,ft::Function,k::Int,q::Int,p::Int)
      new(μv,μFuncR(fr,k,q),μFuncT(ft,p))
   end
   function μOb(μv::Float64,fr::μFuncR,ft::μFuncT)
      new(μv,fr,ft)
   end
end

function eval_μop_r(op::μFuncR,ψb::RPsi,ψk::RPsi)::SparseMatrixCSC{Float64,Int} 
   op.f(ψb,ψk,op.k,op.q)
end
function eval_μop_t(op::μFuncT,ψb::TPsi,ψk::TPsi)::SparseMatrixCSC{Float64,Int} 
   op.f(ψb,ψk,op.k,op.q)
end

function int_enact_1(μf,ψb,ψk)::SparseMatrixCSC{Float64,Int}
   return kron(μf.v,
               eval_μop_t(μf.ft,ψb.T,ψk.T),
               eval_μop_r(μf.fr,ψb.R,ψk.R))
end
function int_enact_2(μf,ψb,ψk,tvecs)::SparseMatrixCSC{Float64,Int}
   return kron(μf.v,
      sparse(tvecs' * eval_μop_t(μf.ft,ψb.T,ψk.T) * tvecs),
         eval_μop_r(μf.fr,ψb.R,ψk.R))
end

function bld_μmat_1s(inthres,μf,ψb,ψk,bvec,kvec)::SparseMatrixCSC{Float64,Int}
   out = int_enact_1(μf[1],ψb,ψk)
   for i ∈ 2:length(μf)
      out += int_enact_1(μf[i],ψb,ψk)
   end
   out = sparse(bvec' * out * kvec)
   droptol!(out,inthres)
   return out
end

function bld_μmat_2s(inthres,μf,ψb,ψk,tvecs,bvec,kvec)::SparseMatrixCSC{Float64,Int}
   out = int_enact_2(μf[1],ψb,ψk,tvecs)
   for i ∈ 2:length(μf)
      out += int_enact_2(μf[i],ψb,ψk,tvecs)
   end
   out = sparse(bvec' * out * kvec)
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
function jvlinds(j,s,vtm)
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return collect(snd:fnd)
end

function intcalc(ctrl,vecs,tvecs,μf,σ)
   s = ctrl.S 
   vtm = ctrl.vtmax 
   mc = ctrl.mcalc 
   nf = ctrl.NFOLD 
   jmax = ctrl.Jmax 
   stages = ctrl.stages 
   mΔj = 1
   jbjk = jbjklister(0.5*iseven(Int(2*s+1)),jmax,mΔj)
   rmsd = Int((2*s+1)*(2*ctrl.mcalc +1))
   cmsd = Int((2*s+1)*(2*ctrl.mcalc +1))
   ints = spzeros(size(vecs,2),size(vecs,2))

   for i in 1:size(jbjk,1)
      jb = jbjk[i,1]
      jk = jbjk[i,2]
      #needs to be corrected for the vt scheme
      kinds = jvlinds(jk,s,vtm)
      binds = jvlinds(jb,s,vtm)
      #filter vecs
      kvecs = vecs[1:Int(2*jk+1)*cmsd,kinds]
      bvecs = vecs[1:Int(2*jb+1)*rmsd,binds]
      ψb = Psi( RPsi(jb,s), TPsi(nf,σ,mc) )
      ψk = Psi( RPsi(jk,s), TPsi(nf,σ,mc) )
      #calculate intensities
      if stages==1
         μs = bld_μmat_1s(ctrl.INTthres ,μf,ψb,ψk,bvecs,kvecs)
      elseif stages==2
         μs = bld_μmat_ws(ctrl.INTthres ,μf,ψb,ψk,bvecs,kvecs,tvecs)
      else
         @warn "this stage count $(stages) not implemented"
      end
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

function freq_gen(min,max,ints,vals,rinds,cinds)
   frqs = spzeros(length(vals),length(vals))
   for i in 1:nnz(ints)
      r = rinds[i]
      c = cinds[i]
      ν = vals[r] - vals[c]
      frqs[r,c] = ν*(min≤abs(ν*0.001)≤max)
   end
   droptol!(frqs,1e-3) #drops all frequencies below 1kHz
   #assign Elow & qns plus sign correction
   rinds, cinds, νs = findnz(frqs)
   return rinds,cinds,νs
end

function qrtcalc(fvals::Array{Float64},TK::Float64)::Float64
   β = -1.0/(20836.61912*TK) #K/MHz
   out = exp.(fvals .* β)
   if size(out,2) > 1
      out[:,2:end] .*= 2.0
   end
   return sum(out)
end

#intensity from B&J: 8π^3 N_A ν * e^(-E_l/kbT) * (1 - e^(-hcν/kbT)) |⟨μ⟩|^2 / 4πϵ0 3hc Q
#Q = ∑_i g_i e^(-E_i / kbT)
#g_i = (2S+1)(2J+1)(1 + σ>0)
function σdegen(σ::Int)::Int
   return !iszero(σ) + 1
end
function σdegen(σ::Vector{Int})::Int
   return prod(σdegen.(σ))
end


function tracalc(ctrl,vals,vecs,tvecs,qns,μf,σ,Qrt)
   ints = intcalc(ctrl,vecs,tvecs,μf,σ)
   rinds, cinds, = findnz(ints)
   rinds, cinds, νs = freq_gen(ctrl.νmin ,ctrl.νmax ,ints,vals,rinds,cinds)
   @show size(νs)
   len = length(νs)
   kbT = ctrl.TK *20836.61912 #MHz/K
   outfrq = zeros(len,4)
   outqns = zeros(Int,len,12)
   degen = (2*ctrl.S +1)*σdegen(σ)
   for i in 1:len
      ν = νs[i]
      r = rinds[i]
      c = cinds[i]
      if ν > 0.0
         outfrq[i,1] = ν
         outfrq[i,3] = vals[c] / csl
         #outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-vals[c]/kbT) - exp(-vals[r]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*degen*(qns[c,1]+1)
         outqns[i,1:6] = qns[r,:]
         outqns[i,7:12] = qns[c,:]
      elseif ν < 0.0
         outfrq[i,1] = -ν
         outfrq[i,3] = vals[r] /csl
         #outfrq[i,4] = √abs(uncs[r,r]^2 + uncs[c,c]^2 - 2*uncs[r,c])
         thermfact = abs(exp(-vals[r]/kbT) - exp(-vals[c]/kbT))/Qrt
         outfrq[i,2] = ints[r,c]*thermfact*degen*(qns[r,1]+1)
         outqns[i,1:6] = qns[c,:]
         outqns[i,7:12] = qns[r,:]
      end
   end
   check = outfrq[:,2] .> ctrl.INTthres 
   outfrq = outfrq[check,:]
   outqns = outqns[check,:]
   #println(outfrq)
   #println(outqns)
   return outfrq, outqns
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