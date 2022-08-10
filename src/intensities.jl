################################################################################
############                  tracalc                   ############
################################################################################


function T(μ,q)
"""
This sloppy function calls the elements of the Cₛ dipole spherical tensor from
   global μa and μb. This will need to be reworked.
"""
   q = Int(q)
   if q==-1
      return (μ[2]-μ[3])/√(2)
   elseif q==0
      return μ[1]
   elseif q==1
      return -(μ[2]+μ[3])/√(2)
   else
      return 0.0
   end
end

function intelem(jb,nb,kb,s,j,n,k,μ)
"""
This implements the dipole moment operator in spherical tensor notation. It is
   an implementation of equation 24-25 of http://dx.doi.org/10.1063/1.1545441
"""
   q = kb-k
   out = wig6j(n,j,s,jb,nb,1)*(-1≤q≤1)*T(μ,q)
   if out == 0.0
      return out
   else
      out = √((2*jb+1)*(2*j+1))*wig3j(n,1,nb,k,q,-kb)*(-1)^(n+s+jb+1)
      out *= √((2*nb+1)*(2*n+1))*(-1)^(n-kb-1)
      return out
   end
end
function intmat(jb,jk,s,mcalc,σ)
"""
This builds the matrix of dipole elements and wang transforms it for the A states
"""
   jdb = Int((2*jb+1)*(2*s+1))
   jdk = Int((2*jk+1)*(2*s+1))
   mat = spzeros(Float64,jdk,jdb)
   omat = spzeros(Float64,jdk,jdb)
   qnb = qngen(jb,s,mcalc,σ)
   qnk = qngen(jk,s,mcalc,σ)
   for y in 1:jdb
   for x in 1:jdk
      @inbounds  mat[x,y] = intelem(jb,qnb[y,2],qnb[y,3],s,jk,qnk[x,2],qnk[x,3],μ[:,1])
      @inbounds omat[x,y] = intelem(jb,qnb[y,2],qnb[y,3],s,jk,qnk[x,2],qnk[x,3],μ[:,2])
   end
   end
   omat = kron(spdiagm(1=>ones(Float64,Int(2*mcalc)),-1=>ones(Float64,Int(2*mcalc))),omat)
   mat = kron(eye(Int(2*mcalc+1)),mat) + omat
   if σ==0||σ==0.0
      Uk = ur(jk,s,mcalc)*ut(mcalc,jk,s)
      Ub = ur(jb,s,mcalc)*ut(mcalc,jb,s)
      mat = Uk*mat*Ub
   else
      Uk = ur(jk,s,mcalc)
      Ub = ur(jb,s,mcalc)
      mat = Uk*mat*Ub
   end
   return mat
end

function intcalc(nmax,s,mcalc,σ,qns,vecs)
   len = size(qns)[1]
   ints = spzeros(Float64,len,len,2)
   for i in 1:len
   jk = qns[i,1]
   lenk = convert(Int,(2*s+1)*(2*jk+1)*(2*mcalc+1))
   for j in i:len
      #fill in int mat
      μmat = intmat(jb,jk,s,mcalc,σ)
      jb = qns[j,1]
      lenb = convert(Int,(2*s+1)*(2*jb+1)*(2*mcalc+1))
      ints[i,j,1] = (transpose(vecs[1:lenk,j])*μmat*vecs[1:lenb,i])^2
   end#for j
   end#for i
   ints[ints .< INTTHRESHOLD] .= 0.0
   ints = dropzeros(ints)
   return ints
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
   trans = zeros(Float64,0,15)
   if isodd(Int64(2*s+1))
      ojb = 1.
      ojk = 0.
   else
      ojb = 1.5
      ojk = 0.5
   end
   lenb = convert(Int,(2*s+1)*(2*ojb+1)*(2*mcalc+1))
   lenk = convert(Int,(2*s+1)*(2*ojk+1)*(2*mcalc+1))
   μmat = intmat(ojb,ojk,s,mcalc,σ)
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
            μmat = intmat(jb,jk,s,mcalc,σ)
         end
         ojb = jb
         ojk = jk
         if (freq > 0.0)&&(freq < 80.0E+03)#&&(Δk < 1)&&(qns[j,3]==0.0)
            int = (transpose(vecs[1:lenk,j])*μmat*vecs[1:lenb,i])^2# *exp(vals[i]/(TK*2.083661912e+4))
            if int>INTTHRESHOLD#&&(abs(qns[j,3])==0.0)&&(abs(qns[i,3])==0.0)
               temp = [freq int vals[i]/csl transpose(qns[j,:]) transpose(qns[i,:])]
               trans = vcat(trans,temp)
            end
         elseif (freq < 0.0)&&(freq > -80.0E+03)#&&(Δk < 1)&&(qns[j,3]==0.0)
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

function tracalc2(nmax,s,mcalc,σ,qns,vals,vecs)
   frqs = spzeros(Float64,length(vals),length(vals))
   ints = spzeros(Float64,length(vals),length(vals))
   MINFREQ = 26500.0
   MAXFREQ = 40000.0
   frqs, ints = intstrcalc!(frqs,ints,mcalc,s,σ,vals,vecs,qns,MAXFREQ,MINFREQ)
   #thermo effects
   #final strength filter
   ints[abs.(ints[:,:]) .< INTTHRESHOLD] .= 0.0
   ints = dropzeros(ints)
   #rearrange
   hind,vind,vs = findnz(ints[:,:])
   out = zeros(Float64,length(hind),16)
   Threads.@threads for i in 1:length(hind)
      #so uhhh I'm not actually sure which is upper state & which is the lower
      h = hind[i]
      v = vind[i]
      out[i,1] = frqs[h,v] #freq
      out[i,2] = ints[h,v] #int
      out[i,3] = vals[h]
      out[i,4] = vals[v]
      out[i,5:end] = [transpose(qns[h,:]) transpose(qns[v,:])]
   end
   return out
end

function intstrcalc!(frqs,ints,mcalc,s,σ,vals,vecs,qns,MAXFREQ,MINFREQ)
   smd = Int((2*s+1)*(2*mcalc+1))
   tind = traindgen(vals)
   Threads.@threads for i in 1:size(tind)[1]
      kind = tind[i,1]
      bind = tind[i,2]
      jb = qns[bind,1]
      jk = qns[kind,1]
      if abs(jb-jk) ≤ 1
         lenb = convert(Int, smd*(2*jb+1))
         lenk = convert(Int, smd*(2*jk+1))
         freq = vals[bind] - vals[kind]
         μmat = intmat(jb,jk,s,mcalc,σ)
         if MAXFREQ ≥ freq ≥ MINFREQ
            frqs[kind,bind] = freq
            ints[kind,bind] = (transpose(vecs[1:lenk,bind])*μmat*vecs[1:lenb,kind])^2
         elseif -MINFREQ ≥ freq ≥ -MAXFREQ
            frqs[bind,kind] = -1.0*freq
            ints[bind,kind] = (transpose(vecs[1:lenk,bind])*μmat*vecs[1:lenb,kind])^2
         else
            #nothing
         end
      end #if
   end #for
   return frqs, ints
end

function lnstrcalc!(trans,smd,σ,vals,qns)

   return trans
end


function traindgen(vals)
   lv = length(vals)
   cv = 1
   inds = zeros(Int,0,2)
   for i in 1:lv
      inds = vcat(inds,hcat(collect(1:cv-1), fill(cv,cv-1)))
      cv += 1
   end
   perm = sortperm(inds[:,1])
   return inds[perm,:]
end
