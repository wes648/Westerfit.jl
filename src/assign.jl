
function assign(nmax,j,s,σ,mcalc,mmax,vals,vecs,rotvs)
   #determine which indices we want to assign
   #nlist = Δlist(j,s)
   #ndegns = @. 2*nlist + 1
   #mcind = mcalc+1
   #mcd = 2*mcalc+1
   #offset = zeros(Int64,length(nlist))
   #offset[2:end] = (mcd) .* ndegns[1:end-1]
   #goals = zeros(Int,0,1)
   #off = 0
   #for n in nlist
   #   start = off + Int((mcalc - mmax)*(2*n+1)) + 1
   #   finsh = off + Int((mcalc + mmax + 1)*(2*n+1))
   #   part = collect(start:finsh)
   #   goals = vcat(goals,part)
   #   off += Int((2*mcalc+1)*(2*n+1))
   #end
   assigned = leadfact(vals,copy(vecs))
   vals = vals[assigned]
   vecs = rotvs*vecs
   vecs = vecs[:,assigned]
   #vals = vals[goals]
   #vecs = vecs[:,goals]
   QNs = qngen(j,s,mcalc,σ)
   return QNs, vals, vecs
end
function leadfact(vals,vecs)
   vecs = abs.(vecs)
   c = zeros(Int,length(vals))#,2)
   for i in 1:length(vals)
      t = argmax(vecs)
      s = t[1] #state number that this eigenvalue maps to
      k = t[2] #eigenvalue number from energetic sorting
      vecs[s,k] = 0.0
      vecs[s,:] = zeros(Float64,size(vecs[1,:])) #prevents state from being double assigned
      vecs[:,k] = zeros(Float64,size(vecs[:,1]))# 3, 4, 1, 2, 5
      c[k] = s
   end
   perm = sortperm(c)
   return perm
end

function suboverlap(vals, tvecs)
   submat = zeros(Float64,(2*mcalc+1),size(tvecs)[2])
   mcounters = zeros(Int,(2*mcalc+1))
   minds = zeros(Int,)
   jd = Int((2*s+1)*(2*j+1))
   si = 1
   fi = jd
   for m in 1:(2*mcalc+1)
      submat[m,:] = sum(tvecs[si:fi,:] .^2, dims=2)
      mcounters[m] = sum(asgn[si:fi] .> 0)
      si += jd
      fi += jd
   end#for
   for i in 1:length(vals)
      t = argmax(submat)
      s = t[1]#m state dominating this eigenvalue
      k = t[2]#eigenvalue number from energetic sorting
      if mcounters[s] < jd
         submat[:,k] = 0.0 #avoid double assignment
#         ind = argmax(tvecs[minds[m],k])
      end#if
   end#for
end#func

#function leadfact(vals,vecs,pr,j,s,m,σ)
#   vecs = abs.(vecs)
#   c = zeros(Int,length(vals))#,2)
#   qns = qngen(j,s,m,σ)
#   for i in 1:length(vals)
#      tvec = vecs[:,i]
#      unresolved=true
#      while unresolved
#         t = argmax(tvec)#state number that this eigenvalue maps to
#         echeck = Ediag(pr,j,s,qns[t,2],qns[t,3],qns[t,4],σ)
#         if (t∉c)#tolcomp(vals[i],echeck)&&
#            vecs[t,:] = zeros(Float64,size(vecs[1,:])) #prevents state from being double assigned
#            vecs[:,i] = zeros(Float64,size(vecs[:,1]))# 3, 4, 1, 2, 5
#            c[i] = t
#            unresolved=false
#         elseif tvec == zeros(size(tvec))
#            println("FUCK j=$j i=$i")
#            break
#         else
#            tvec[t] = 0.0
#         end#if
#      end#while
#   end#for
#   perm = sortperm(c)
#   return perm
#end

function tolcomp(a,b)
   #is a ∈ b±0.5tol
   tol = 0.10
   lw = b*(1 - 0.5*tol)
   up = b*(1 + 0.5*tol)
   test = (a ≥ lw)&&(a ≤ up)
end

function assign3(nmax,j,s,σ,mcalc,mmax,vals,vecs)
#This is the 3 layer qn assignment routine
#define goal states
#leadfact layer
#suboverlap layer
#ondiag-comp layer
end

function inchecker(a,b)
   test = true
   for i in 1:length(a)
      test *= a[i] ∈ b
   end
   return test
end

function leadfactlayer(vals,vecs)
   tvecs = abs.(copy(vecs))
   c = zeros(Int,length(vals))#,2)
   unassigned = length(vals)
   for i in 1:length(vals)
      t = argmax(tvecs)
      elem = maximum(tvecs)
      if elem > 0.7
         s = t[1] #state number that this eigenvalue maps to
         k = t[2] #eigenvalue number from energetic sorting
         tvecs[s,:] = zeros(Float64,size(tvecs[1,:])) #prevents state from being double assigned
         tvecs[:,k] = zeros(Float64,size(tvecs[:,2]))# 3, 4, 1, 2, 5
         c[k] = s
         unassigned -= 1
      else
         println("Moving to suboverlap layer")
         break
      end#if
   end#for
   #perm = sortperm(c)
   return c, tvecs#perm
end#func

function suboverlaplayer(asgn, vals, vecs)
   unassigned = unique([collect(1:length(vals)); asgn])
   #this uses the filtered eigvenvectors from the previous layer
   tvecs = abs.(vecs) .^2
   submat = zeros(Float64,(2*mcalc+1),size(tvecs)[2])
   mcounters = zeros(Int,(2*mcalc+1))
   minds = zeros(Int,)
   jd = Int((2*s+1)*(2*j+1))
   si = 1
   fi = jd
   for m in 1:(2*mcalc+1)
      submat[m,:] = sum(tvecs[si:fi,:].^2, dims=2)
      mcounters[m] = sum(asgn[si:fi] .> 0)
      si += jd
      fi += jd
   end#for
   for i in 1:length(vals)
      t = argmax(submat)
      s = t[1]#m state dominating this eigenvalue
      k = t[2]#eigenvalue number from energetic sorting
      if mcounters[s] < jd
         submat[s,k] = 0.0
#         ind = argmax(tvecs[minds[m],k])
      end#if
   end#for
end#func
