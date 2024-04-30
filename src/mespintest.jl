using LinearAlgebra, SparseArrays, WIGXJPFjl
#each A block has 8 states, each E 4
qns = [ 
 0  1.5  1  -1  -3  0 #1
 0  1.5  1   0  -3  0
 0  1.5  1   1  -3  0
 0  1.5  2  -2  -3  0
 0  1.5  2  -1  -3  0
 0  1.5  2   0  -3  0
 0  1.5  2   1  -3  0
 0  1.5  2   2  -3  0 #8
 0  0.5  0   0  -2  1 #9
 0  0.5  1  -1  -2  1 
 0  0.5  1   0  -2  1 
 0  0.5  1   1  -2  1 #12
 0  0.5  0   0  -1  2 #13
 0  0.5  1  -1  -1  2
 0  0.5  1   0  -1  2
 0  0.5  1   1  -1  2 #16
 0  1.5  1  -1   0  0 #17
 0  1.5  1   0   0  0
 0  1.5  1   1   0  0
 0  1.5  2  -2   0  0
 0  1.5  2  -1   0  0
 0  1.5  2   0   0  0
 0  1.5  2   1   0  0
 0  1.5  2   2   0  0 #24
 0  0.5  0   0   1  1 #25
 0  0.5  1  -1   1  1 
 0  0.5  1   0   1  1
 0  0.5  1   1   1  1 #28
 0  0.5  0   0   2  2 #29
 0  0.5  1  -1   2  2 
 0  0.5  1   0   2  2
 0  0.5  1   1   2  2 #32
 0  1.5  1  -1   3  0 #33
 0  1.5  1   0   3  0
 0  1.5  1   1   3  0 
 0  1.5  2  -2   3  0
 0  1.5  2  -1   3  0
 0  1.5  2   0   3  0
 0  1.5  2   1   3  0
 0  1.5  2   2   3  0] #40
#F   J   N   K   m  Γ
#1   2   3   4   5  6
function ur(n::Int)::SparseMatrixCSC{Float64, Int64}
      c = spzeros(n)
      r = c'
      u = sparse((1/√2)*I,n,n)
      l = rotl90(u)
      out = [-u c l; r 1 r; l c u]
   return out
end
function utbuild()
   out = spzeros(40,40)
   out += I(40)
   out[1:16,1:16] *= -1
   out[1:3,33:35] = I(3) #±3 link
   out[4:8,36:40] = I(5)
   out[33:35,1:3] = I(3)
   out[36:40,4:8] = I(5)
   out[9,25] = 1 #±1 link
   out[25,9] = 1 
   out[10:12,26:28] = I(3)
   out[26:28,10:12] = I(3) 
   out[13,29] = 1 #±2 link
   out[29,13] = 1
   out[14:16,30:32] = I(3)
   out[30:32,14:16] = I(3)
   out *= 1/√2
   out[17:24,17:24] = I(8)
   return sparse(out)
end
function urbuild()
   out = ur(1)
   nlist = [1;2;0;1;0;1;1;2;0;1;0;1;1;2]
   for n in nlist[2:end]
      out = cat(out,ur(n),dims=(1,2))
   end
   return sparse(out)
end

δ(x::Number,y::Number)::Float64 = x==y
□rt(x::Number)::Float64 = √(x*(x>zero(x)))
fh(x::Number,y::Number)::Float64 = □rt((x-y)*(x+y+1))

function rotop(k,q,nb,kb,nk,kk)
   out = δ(nb,nk)*(-1)^(nb-kb)
   out *= wig3j(nb,2,nk,-kb,q,kk)
   out *= nk*(nk+1)*(2*nk+1)*wig6j(nk,nk,1,k,1,nb)
   out *= √(2*k+1)*(-1)^k
   return out
end
function hrot(qnb,qnk)
   T00 = 1.0
   T2r = [1.0; -1.0; 1.0; 1.0; 1.0]
   jb = qnb[2]
   jk = qnk[2]
   mb = qnb[5]
   mk = qnk[5]
   if (jb==jk)&&(mb==mk)
      nb = qnb[3]
      kb = qnb[4]
      nk = qnk[3]
      kk = qnk[4]
      out = rotop(0,0,nb,kb,nk,kk)*T00
      for q in -2:2
         out += rotop(2,q,nb,kb,nk,kk)*T2r[q+3]
      end
   else
      out = 0.0
   end
   return out
end

function srop(k,q,j,s,nb,kb,nk,kk)
   out  = √(nk*(nk+1)*(2*nk+1))*wig6j(1,1,k,nb,nk,nk)*(-1)^k
   out += √(nb*(nb+1)*(2*nb+1))*wig6j(1,1,k,nk,nb,nb)
   out *= √((2*k+1)*s*(s+1)*(2*s+1)*(2*nk+1)*(2*nb+1))*(-1)^(j+s+nb)
   out *= wig6j(nk,s,j,s,nb,1)*0.5
   out *= wig3j(nb,k,nk,-kb,q,kk)*(-1)^(nb-kb)
   return out
end
function hsr(qnb,qnk)
   s = 1/2
   T00 = 1.0
   T1sr = [1.0; 0.0; 1.0]
   T2sr = [1.0; -1.0; 1.0; 1.0; 1.0]
   jb = qnb[2]
   jk = qnk[2]
   mb = qnb[5]
   mk = qnk[5]
   if (jb==jk)&&(mb==mk)
      nb = qnb[3]
      kb = qnb[4]
      nk = qnk[3]
      kk = qnk[4]
      out = srop(0,0,jk,s,nb,kb,nk,kk)*T00
      for q in -1:1
         out += srop(1,q,jk,s,nb,kb,nk,kk)*T1sr[q+2]
      end
      for q in -2:2
         out += srop(2,q,jk,s,nb,kb,nk,kk)*T2sr[q+3]
      end
   else
      out = 0.0
   end
   return out
end

function quaop(q,f,i,s,jb,nb,kb,jk,nk,kk)
   out = 1/wig3j(i,2,i,-i,0,i)
   out *= (-1)^(jk+i+f+nb+s+jb) *wig6j(i,jb,f,jk,i,2)*wig6j(nb,jb,s,jk,nk,2)
   out *= √((2*jb+1)*(2*jk+1)*(2*nb+1)*(2*nk+1))*(-1)^(nb-kb)
   out *= wig3j(nb,2,nk,-kb,q,kk)
   return out
end
function hqua(qnb,qnk)
   s = 1/2
   i = 3/2
   T2q = [1.0; -1.0; 1.0; 1.0; 1.0]
   mb = qnb[5]
   mk = qnk[5]
   Γk = qnk[6]
   out = 0.0
   if (mb==mk)&&(Γk==zero(Γk))
      f = qnk[1]
      jb = qnb[2]
      nb = qnb[3]
      kb = qnb[4]
      jk = qnk[2]
      nk = qnk[3]
      kk = qnk[4]
      for q in -2:2
         out += quaop(q,f,i,s,jb,nb,kb,jk,nk,kk)*T2q[q+3]
      end
   end
   return out
end

function ired(ib,Γb,ik,Γk)
   if Γb==Γk
      out = √(ik*(ik+1)*(2*ik+1))
   elseif (Γb+Γk)==3
      out = -2*√(ik*(ik+1)*(2*ik+1))
   else
      out = √6
   end
   return out
end
function fermiop(f,s,n,jb,Γb,jk,Γk)
   ib = 1/2 + (Γb==zero(Γb))
   ik = 1/2 + (Γk==zero(Γk))
   out = (-1)^(n+s+jb+jk+ib+f+1)
   out *= √((2jb+1)*(2jk+1)*s*(s+1)*(2s+1))*wig6j(ib,jb,f,jk,ik,1)
   out *= wig6j(s,jb,n,jk,s,1)*ired(ib,Γb,ik,Γk)
   return out
end
function spspop(q,f,s,jb,nb,kb,Γb,jk,nk,kk,Γk)
   ib = 1/2 + (Γb==zero(Γb))
   ik = 1/2 + (Γk==zero(Γk))
   out = √(30*s*(s+1)*(2s+1)*(2jb+1)*(2jk+1)*(2nb+1)*(2nk+1) )*(-1)^(jk+ib+f) 
   out *= wig6j(ib,jb,f,jk,ik,1)*wig9j(nb,nk,2,s,s,1,jb,jk,1)*ired(ib,Γb,ik,Γk)
   out *= wig3j(nb,2,nk,-kb,q,kk)*(-1)^(nb-kb)
   return out
end
function hspsp(qnb,qnk)
   s = 1/2
   af = 1.0
   T2s = [1.0; -1.0; 1.0; 1.0; 1.0]
   af = 1.0
   T2s = zeros(5)
   mb = qnb[5]
   mk = qnk[5]
   out = 0.0   
   if (mb==mk)
      f = qnk[1]
      jb = qnb[2]
      nb = qnb[3]
      kb = qnb[4]
      Γb = qnb[6]
      jk = qnk[2]
      nk = qnk[3]
      kk = qnk[4]
      Γk = qnk[6]
      if (δ(nb,nk)==1)&&((δ(kb,kk)==1))
         out += fermiop(f,s,nk,jb,Γb,jk,Γk)*af
      end
      for q in -2:2
         #out += spspop(q,f,s,jb,nb,kb,Γb,jk,nk,kk,Γk)*T2s[q+3]
      end
   end
   return out
end

function htor(qnb,qnk)
   F = 1.0
   V = 1.0
   ρ = 1.0
   jb = qnb[2]
   nb = qnb[3]
   kb = qnb[4]
   mb = qnb[5]
   jk = qnk[2]
   nk = qnk[3]
   kk = qnk[4]
   mk = qnk[5]
   out = δ(jb,jk)*δ(nb,nk)*(δ(kb,kk)*(
      δ(mb,mk)*(F*mk^2 - 2F*ρ*mk*kk + V/2) +
      δ(mb,mk+3)*(-V/4) + δ(mb,mk-3)*(-V/4) )
      + δ(mb,mk)*δ(kb,kk+1)*mk*fh(nk,kk) 
      + δ(mb,mk)*δ(kb,kk-1)*mk*fh(nk,kb)
      )
   return out
end

function hbuild(qns)
   out = spzeros(size(qns,1),size(qns,1))
   for a in 1:size(qns,1)
      qnk = qns[a,:]
   for b in 1:size(qns,1)
      qnb = qns[b,:]
      out[a,b] += hrot(qnb,qnk)
      out[a,b] += hsr(qnb,qnk)
      out[a,b] += hqua(qnb,qnk)
      out[a,b] += hspsp(qnb,qnk)
      out[a,b] += htor(qnb,qnk)      
   end
   end
   return out
end

U = urbuild()*utbuild()
H = hbuild(qns)