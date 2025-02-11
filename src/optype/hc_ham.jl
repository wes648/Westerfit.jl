
hr2on(ns,ks,bk::Float64,bn::Float64) = @. bn*eh2(ns) + bk*ks^2 
hr2of1(ns,ks,dab::Float64) = @. dab*(ks-0.5)*fh(ns,ks-1)
hr2of2(ns,ks,bpm::Float64) = @. bpm*fh(ns,ks-1)*fh(ns,ks-2)
function hrot2(pr::Vector{Float64},ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   out = spzeros(size(ns,1),size(ns,1))
   #p0 = hr2on(ns,ks,pr[1],pr[2])
   out[diagind(out)] .= hr2on(ψ.N,ψ.K,pr[1],pr[2])
   #p1 = hr2of1(ns[2:end],ks[2:end], pr[4])
   out[diagind(out,1)] .= hr2of1(ψ.N[2:end],ψ.K[2:end], pr[4])
   #p2 = hr2of2(ns[3:end],ks[3:end], pr[3])
   out[diagind(out,2)] .= hr2of2(ψ.N[3:end],ψ.K[3:end], pr[3])
   #out = spdiagm(0=>p0,1=>p1,2=>p2)
   return dropzeros!(out)
end


function nsred(l::Int,nb::Int,nk::Int)::Float64
   if (abs(l)==1)&&(nb==nk)
      out = 0.0
   else
   out = 0.5*( 
   √(nk*(nk+1)*(2*nk+1))*
      wig6j( 1, 1, l,
            nb,nk,nk) + #*powneg1(l) + 
   √(nb*(nb+1)*(2*nb+1))*
      wig6j( 1, 1, l,
            nk,nb,nb))
   end
   return out
end
function jsred(j,s,nb::Int,nk::Int)::Float64
   return wig6j(nk, s, j,
                 s,nb, 1)*√((2*nb+1)*(2*nk+1))
end
function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)::Float64
   return pr*wig3j( nb,l,nk,
                   -kb,q,kk)*√(2*l+1)*
       nsred(l,nb,nk)*jsred(j,s,nb,nk)*powneg1(nb-nk-kb)
end
function hsr2(pr::Array{Float64},j,s,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   le = convert(Int,(2j+1)*(2s+1))
   out = spzeros(le,le)
   ns = srinds(j,s)
   for i in 1:length(ns), j in i:length(ns)
      nb = ns[j]
      nk = ns[i]
      nib = ni[j] + nb
      nik = ni[i] + nk
      fact = nsred(l,nb,nk)*jsred(j,s,nb,nk)*√(2l+1)*powneg1(nb-nk)
      for a in -nk:nk, b in -nb:nb
         q = a - b
         if abs(q) < 3
out[nib+b,nik+a] = pr*fact*wig3j(nb,l,nk,
                                 -b,q, a)*powneg1(kb)
         end#if 
      end#for ks
   end#for ns
end
function hsr(pr::Array{Float64},j,s,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   le = length(ψ.K)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   ts = SA[0  1 1 2  2 2  2 2;
           0 -1 1 0 -1 1 -2 2;
           1  2 2 3  4 4  5 5;
           1  1 1 1 -1 1  1 1]
   for i in 1:8
      tv = ts[:,i]
      prm = pr[tv[3]]*tv[4]
      if prm ≠ 0.0
      for a in 1:le, b in a:le
         nb = ψ.N[b]
         kb = ψ.K[b]
         nk = ψ.N[a]
         kk = ψ.K[a]
         if abs(nb-nk)≤1 && (tv[2]+kk-kb)==0
            out[b,a] += srelem(prm,tv[1],tv[2], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#sr ind for loop
      end#prm chck if
   end#sr term for loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(j+s)
   return out
end

function wiginv(s::Real)::Float64
   return s<one(s) ? 0.0 : inv(wig3j( s,2,s,
                                     -s,0,s))
end
function qured(j,s,nb,nk)
   return jnred(nb,nk)*wig6j(j, s,nb,
                             2,nk, s)
end
function qulm(pr,q,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   return pr*qured(j,s,nb,nk)*
#             δ(nb,nk)* #This line can be used to emulate the perturbative 
             wig3j( nb, 2,nk,
                   -kb, q,kk)*powneg1(nb+nk-kb)
end
#const tq = SA[0 -1 1 -2 2; 
#              1  2 2  3 3;
#              1 -1 1  1 1]
function hqu(pr,j,s,ψ::Psi)::SparseMatrixCSC{Float64, Int64}
   le = size(ψ.K,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   tq = SA[0 -1 1 -2 2; 
           1  2 2  3 3;
           1 -1 1  1 1]
   for i in 1:5
      tv = tq[:,i]
      prm = pr[tv[2]]*tv[3]
      if prm ≠ 0.0
      for a in 1:le, b in a:le
         nb = ψ.N[b]
         kb = ψ.K[b]
         nk = ψ.N[a]
         kk = ψ.K[a]
         if abs(nb-nk)≤2 && (tv[1]+kk-kb)==0
            out[b,a] += qulm(prm,tv[1], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#qu ind for loop
      end#prm chck if
   end#qu term for loop
   dropzeros!(out)
   out .*= 0.25*wiginv(s)*powneg1(j+s)
   #@show out
   return out
end
