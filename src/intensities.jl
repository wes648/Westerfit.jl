################################################################################
############                  tracalc                   ############
################################################################################


function Tμ(q)
"""
This sloppy function calls the elements of the Cₛ dipole spherical tensor from
   global μa and μb. This will need to be reworked.
"""
   q = Int(q)
   if q==-1
      return μb/√(2)
   elseif q==0
      return μa
   elseif q==1
      return -μb/√(2)
   else
      return 0.0
   end
end

function intelem(jb,nb,kb,s,j,n,k)
"""
This implements the dipole moment operator in spherical tensor notation. It is
   an implementation of equation 24-25 of http://dx.doi.org/10.1063/1.1545441
"""
   q = kb-k
   out = WIGXJPF.wig6j(n,j,s,jb,nb,1)*WIGXJPF.wig3j(n,1,nb,k,q,-kb)*Tμ(q)
   if out == 0.0
      return out
   else
      out = √((2*jb+1)*(2*j+1))*(-1)^(n+s+jb+1)
      out *= √((2*nb+1)*(2*n+1))*(-1)^(n-kb-1)
      return out
   end
end
function intmat(jb,jk,s,mcalc,σ,lenb,lenk)
"""
This builds the matrix of dipole elements and wang transforms it for the A states
"""
   mat = zeros(Float64,lenk,lenb)
   qnb = qngen(jb,s,mcalc,σ)
   qnk = qngen(jk,s,mcalc,σ)
   for y in 1:lenb
   for x in 1:lenk
      @inbounds mat[x,y] = intelem(jb,qnb[y,2],qnb[y,3],s,jk,qnk[x,2],qnb[x,3])
   end
   end
   if σ==0
      Uk = ur(jk,s,mcalc)*ut(mcalc,jk,s)
      Ub = ur(jb,s,mcalc)*ut(mcalc,jb,s)
      mat = Uk*mat*Ub
   end
   return mat
end

function tracalc(nmax,s,mcalc,σ,qns,vals,vecs)
"""
This repulsively slow function calculates all of the transitions from the
   eigenvalues & vectors from tsrdiag. This does include an intensity cutoff but
   the only hard coded selction rule is that |ΔJ| ≤ 1. Every pair that meets that
   term is tested since the meaning of the other quantum numbers is incredibly
   limited.
"""
   jmax = nmax - s
   ns = Δlist(jmax,s)
   trans = zeros(Float64,0,15)
   karray = kron(ones(2*mcalc+1),collect(-ns[1]:ns[1]))
   for i in 2:length(ns)
      karray = vcat(karray,kron(ones(2*mcalc+1),collect(-ns[i]:ns[i])))
   end
   if isodd(Int64(2*s+1))
      ojb = 1.
      ojk = 0.
   else
      ojb = 1.5
      ojk = 0.5
   end
   lenb = convert(Int,(2*s+1)*(2*ojb+1)*(2*mcalc+1))
   lenk = convert(Int,(2*s+1)*(2*ojk+1)*(2*mcalc+1))
   μmat = intmat(ojb,ojk,s,mcalc,σ,lenb,lenk)
   for i in 1:length(vals)
   for j in (i+1):length(vals)
      Δj = abs(qns[j,1] - qns[i,1])
      Δk = abs(qns[j,3] - qns[i,3])
      if Δj ≤ 2
#      Δka = abs(qns[j,3] - qns[i,3])
#      if Δka ≤ 2
         jb = qns[j,1]
         jk = qns[i,1]
         lenb = convert(Int,(2*s+1)*(2*jb+1)*(2*mcalc+1))
         lenk = convert(Int,(2*s+1)*(2*jk+1)*(2*mcalc+1))
         freq = vals[j] - vals[i]
         if (jb!=ojb)||(jk!=ojk)
            μmat = intmat(jb,jk,s,mcalc,σ,lenb,lenk)
         end
         ojb = jb
         ojk = jk
         if (freq > 0.0)&&(freq < 80.0E+03)&&(Δk < 1)&&(qns[j,3]==0.0)
            int = (transpose(vecs[1:lenk,j])*μmat*vecs[1:lenb,i])^2# *exp(vals[i]/(TK*2.083661912e+4))
            if int>INTTHRESHOLD#&&(abs(qns[j,3])==0.0)&&(abs(qns[i,3])==0.0)
               temp = [freq int vals[i]/csl transpose(qns[j,:]) transpose(qns[i,:])]
               trans = vcat(trans,temp)
            end
         elseif (freq < 0.0)&&(freq > -80.0E+03)&&(Δk < 1)&&(qns[j,3]==0.0)
            int = (transpose(vecs[1:lenk,i])*μmat*vecs[1:lenb,j])^2# *exp(vals[j]/(TK*2.083661912e+4))
            if (int>INTTHRESHOLD)#&&(abs(qns[j,3])==0.0)&&(abs(qns[i,3])==0.0)
               temp = [-freq int vals[j]/csl transpose(qns[i,:]) transpose(qns[j,:])]
               trans = vcat(trans,temp)
            end
         else
         end#freq if
#      else
#      end#Δka if
      else
         break
      end#Δj if
   end#i for
   end#j for
   return trans
end

function pred2lne(sim)
"""
Converts the transitions output from westersim into the line format for
   westerfit. Allows for quick test scripting
"""
   out = zeros(Float64,size(sim)[1],14)
   out[:,1:12] = sim[:,4:end]
   out[:,13] = sim[:,1]
   out[:,14] = fill(0.08,size(sim)[1])
   return out
end
