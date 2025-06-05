
function enact_init(O::Op,ψ::RPsi,val)::SparseMatrixCSC{Float64,Int}
   #out = O.a≠0 ? O.v .* eh.(ψ.N).^O.a : fill(O.v,ψ.lng)
   out = O.a≠0 ? n2(0.5*val, ψ.N, O.a) : fill(0.5*val,ψ.lng)
   if O.b≠0; out .*= ns(ψ.J,ψ.S,ψ.N,O.b); end
   if O.c≠0; out .*= eh(ψ.S)^O.c ; end
   if O.d≠0; out .*= nz(ψ.K, O.d) ; end
   return spdiagm(out)
end

function torop(a::Int,b::Int,nf::Int,ms::StepRange{Int,Int})::SparseMatrixCSC{Float64,Int}
#performance of this function is deeply lacking
   if !iszero(b)
      @inbounds part = 0.5.* (ms[1+b:end] .- nf*b).^a
      out = spdiagm(b=>part)
      out[diagind(out,-b)] .= reverse!(part * powneg1(a))
      dropzeros!(out)
   elseif !iszero(a) && iszero(b)
      out = spdiagm(0=>ms .^a)
   else
      out = spdiagm(ones(size(ms)))
   end
   #@show out
   return out
end
function enact_tor(tp::Array{Int,2},ψ::TPsi,
                   ut::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   out = torop(tp[1,1],tp[2,1],ψ.nf[1],ψ.ms[1])
   if iszero(ψ.σ[1])
      out = dropzeros!(ut*out*ut)
   end
   @inbounds for i in 2:size(tp,2)
      part = torop(tp[1,i],tp[2,i],ψ.nf[i],ψ.ms[i])
      if iszero(ψ.σ[i])
         part = dropzeros!(ut*out*ut)
      end
      out = kron(part,out)
   end
   return out
end

function enact(O::Op,ψ::Psi,val::Float64,ur::SparseMatrixCSC{Float64,Int},
               ut::SparseMatrixCSC{Float64,Int})::SparseMatrixCSC{Float64,Int}
   out = enact_init(O,ψ.R,val)
   @inbounds for i in 1:length(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   out = dropzeros!(ur*out*ur)
   if !iszero(O.tp) 
      part = enact_tor(O.tp,ψ.T,ut)
      out = kron(part,out)
   else
      out = kron(I(ψ.T.lng),out)
   end
   return out
end
#This allows the basis set to be distributed among a list of added Operators
function enact(O::Vector{Op},ψ::Psi,val::Vector{Float64},stgs::Vector{Int},
         ur::SparseMatrixCSC{Float64,Int},ut::SparseMatrixCSC{Float64,Int}
         )::SparseMatrixCSC{Float64,Int}
   out = enact(O[1],ψ,val[1],ur,ut)
   @inbounds for i in 2:length(O)
      if stgs[i] ≥ 0
         out += enact(O[i],ψ,val[i],ur,ut)
      else
         out += enact(O[i],ψ,val[i],ur,ut)*val[i+stgs[i]]
   end;end
   tplus!(out)
   return out
end

#tsr_diag variants: 0 generic universal, 1 preloads Hrs, 2 two-stage preloads Hrs
function tsrdiag_0(ℋ::Vector{Op},ψ::Psi)
   H = enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
#having U at the end of kron maintians the tight block structure of the lower index tops
   if ψ.σ[i] == 0
      U = kron(ur(mcalc),U)
   else
      U = kron(sparse(1.0I,2mcalc+1,2mcalc+1),U)
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = droptol!(U*H*U,2*eps())
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl.assign =="ram36")||(ctrl.assign =="RAM36")
      perm = ramassign(vecs,j,s,mcalc,vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl.assign =="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl.assign =="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   return vals, vecs
end
function tsrdiag_1(Hr::SparseMatrixCSC{Float64,Int},vals::Vector{Float64},ctrl,stgs,
                  ℋ::Vector{Op},ψ::Psi,σs)
   H = kron(I(ψ.T.lng),Hr) #<--------- CHECK
   H += enact(ℋ,ψ,vals[12:end],stgs,ur(ψ.R.J,ψ.R.S),ur(ctrl.mcalc )) #<--------- CHECK
   #tplus!(H)
   vals,vecs = eigen!(Symmetric(Matrix(H),:L))
   if (ctrl.assign =="ram36")||(ctrl.assign =="RAM36")
      perm = ramassign(vecs,ψ.R.J,ψ.R.S,ctrl.mcalc ,ctrl.vtmax ) #<--------- CHECK
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl.assign =="expectk"
      vals, vecs = expectkassign!(vals,vecs,ψ.J,ψ.S,nf,ctrl.mcalc ,σ)      
   elseif ctrl.assign =="eeo"
      vals, vecs = eeoassign!(vals,vecs,ψ.J,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,ψ.J,s,nf,mcalc,σ)
   end
   return vals, vecs
end
function tsrcalc_1stg!(vals,vecs,jlist,σs,ctrl,prm,stg,ℋ)
   σcnt = size(σs,2)
   msd = (2*ctrl.mcalc +1)*length(ctrl.NFOLD )
   for ind ∈ axes(jlist,1)
      jd = Int(jlist[ind,1]+1)
      sc = jlist[ind,2]
      dest = jvdest2(0.5*jlist[ind,1],ctrl.S ,ctrl.vtmax )
      ϕ = RPsi(0.5*jlist[ind,1],ctrl.S )
      Hrot = hrot2(prm[1:4],ϕ)
      if ctrl.S ≥1.0
         Hrot += hsr(prm[5:8],ψ.J,ψ.S,ϕ) + hqu(prm[9:11],ψ.J,ψ.S,ϕ)
      elseif ctrl.S ==0.5
         Hrot += hsr(prm[5:8],ψ.J,ψ.S,ϕ)
      end
      ψ = Psi(ϕ,TPsi(ctrl.NFOLD ,σs[:,sc],ctrl.mcalc ))
      vals[dest,sc],vecs[1:jd*msd,dest,sc] = tsrdiag_1(Hrot,prm,ctrl,stg,ℋ,ψ,σs[:,sc])
   end#j
   return vals, vecs
end#f

function torcalc!(tvals,tvecs,ctrl,prm,ℋ,ϕ,stg,sc)
   tsize = (2*ctrl.mcalc +1)^length(ctrl.NFOLD )
   dest = 1:ctrl.mmax +1
   H = torbuild(prm[12:end],ℋ,ϕ,stg,tsize,ctrl.mcalc )
   H = Matrix(H)
   H = eigen!(H)
   tvals[:,sc] = H.values[dest]
   tvecs[:,:,sc] = H.vectors[:,dest]
end
function torbuild(vals,O::Vector{Op},ψ::TPsi,stgs,msz,mc)::SparseMatrixCSC{Float64,Int}
   out = spzeros(msz,msz)
   U = ur(mc)
   @inbounds for i in 1:length(O)
      if isone(stgs[i])
         out += enact_tor(O[i].tp,ψ,U)*vals[i]
      elseif (stgs[i] < 0) && isone(stgs[i+stgs[i]])
         out += enact_tor(O[i].tp,ψ,U)*vals[i+stgs[i]]*vals[i]
      else
      end
   end
   #@show out
   return dropzeros!(out)
end
function enact_stg2(O::Op,ψ::Psi,val,tvcs,mc,ur,ut)::SparseMatrixCSC{Float64,Int}
   #printstyled("start\n",color=:green)
   #@show O.v
   out = enact_init(O,ψ.R,val)
   @inbounds for i in 1:length(O.rf)
      out *= eval_rop(O.rf[i],ψ.R)
   end
   #@show out
   if !iszero(O.tp) #O.tp ≠ zeros(Int,size(O.tp))
      part = enact_tor(O.tp,ψ.T,ut)
      part = dropzeros!(sparse(tvcs' * part * tvcs))
      #@show part
      out = kron(part,out)
   else
      out = kron(I(size(tvcs,2)),out)
   end
   if !isdiag(out) 
      tplus!(0.5*out) 
   end
   #@show out
   #printstyled("stop\n",color=:red)
   return out #<- dispatch 
end
function h_stg2build!(Hmat,O::Vector{Op},ψ::Psi,vals,stgs,siz,tvcs,mc
      )::SparseMatrixCSC{Float64,Int}
   Ur = ur(ψ.R.J,ψ.R.S)
   Ut = ur(mc)
   part = spzeros(size(Hmat))
      @inbounds for i in 1:length(O)
      if iszero(stgs[i]) #this is for future oddities 
         part .+= enact_stg2(O[i],ψ,vals[i],tvcs,mc,Ur,Ut)
      elseif stgs[i] < 0 && iszero(stgs[i+stgs[i]])
         part .+= enact_stg2(O[i],ψ,vals[i]*vals[i+stgs[i]],tvcs,mc,Ur,Ut)
      #else
      end
   end
   Hmat .+= tplus!(part)
   return Hmat
end

function tsrdiag_2(Hr::SparseMatrixCSC{Float64,Int},ctrl,tvals,tvecs,ℋ::Vector{Op},
                  ψ::Psi,prm,stg)
   #printstyled("ψ.J = $(ψ.J), ψ.σ = $(ψ.σ)\n",color=:cyan)
   H = kron(I(ctrl.mmax +1)^length(ctrl.NFOLD ),Hr) 
   H[diagind(H)] .+= kron(tvals, ones(Int((2ψ.R.J+1)*(2ψ.R.S+1)) ))
   h_stg2build!(H,ℋ,ψ,prm[12:end],stg,(2*ctrl.mcalc +1)^length(ctrl.NFOLD ),
               tvecs,ctrl.mcalc )
   #tplus!(H)
   #wangtrans2!(H,ctrl.mmax ,ψ)
   vals,vecs = eigen!(Symmetric(Matrix(H),:L))
   #@show vecs
   perm = twostg_assign(vecs,ψ.R.J,ψ.R.S,ctrl.mmax ,ctrl.vtmax )
   vals = vals[perm]
   vecs = vecs[:,perm]
   #@show nnz(dropzeros(sparse(vals)))
   return vals, vecs
end

function tsrcalc_2stg!(vals,vecs,tvals,tvecs,jlist,σs,ctrl,prm,stg,ℋ)
   σcnt = size(σs,2)#this doesn't work for 1 top
   tsize = (2*ctrl.mcalc +1)^length(ctrl.NFOLD )
   for sc in 1:σcnt
      ϕ = TPsi(ctrl.NFOLD ,σs[sc],ctrl.mcalc )
#      @show σs[sc]
      torcalc!(tvals,tvecs,ctrl,prm,ℋ,ϕ,stg,σs[:,sc])
      #tvals,tvecs = eigen!(Matrix(torbuild(ℋ,ϕ,stg,tsize)))
   end
   msd = Int(2ctrl.S +1)*(ctrl.mmax +1)
for j in jlist
   jd = Int(2j+1)
   dest = jvdest2(j,ctrl.S ,ctrl.vtmax ) 
   ψ = RPsi(j,ctrl.S )
   Hrot = hrot2(prm[1:4],ψ)
   if ctrl.S ≥1.0
      Hrot += hqu(prm[9:11],ψ.J,ψ.S,ψ)
      if norm(prm[5:8]) > 0.0
         Hrot += hsr(prm[5:8],ψ.J,ψ.S,ψ)
      end
   elseif ctrl.S ==0.5
      Hrot += hsr(prm[5:8],ψ.J,ψ.S,ψ)
   end
   Hrot = sparse(Symmetric(Hrot,:L))
   for sc in 1:σcnt
      ϕ = Psi(ψ,TPsi(ctrl.NFOLD ,σs[:,sc],ctrl.mcalc ))
      vals[dest,sc],vecs[1:jd*msd,dest,sc] = tsrdiag_2(Hrot,ctrl,tvals[:,sc],tvecs[:,:,sc],
                                                ℋ,ϕ,prm,stg)
   end#σs
end#j
   return vals,vecs,tvals,tvecs
end#f
