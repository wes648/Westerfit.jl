

function englin(s,eng,qunl,σ)
   if iszero(s)
      part = lpad(qunl[2],4)*"," # N
   else
     # part = lpad(qunl[1],4)*"/2," #J
      part = lpad(qunl[1],4)*","
      part *= lpad(qunl[2],4)*"," # N
   end
   part *= lpad(qunl[3],4)*"," # Ka
   part *= lpad(qunl[4],4)*"," # Kc
   part *= lpad(qunl[5],4)*"," # vt
   part *= lpad(    σ-1,4)*"," # σ
   part *= " "*lpad(@sprintf("%0.10f", eng), 16) # eng
   return part
end

function engwriter(molnam, jmax,s,vtm,vals)
   σcnt = size(vals,2)
   io = open("$molnam.eng", "w") do io
      if iszero(s)
      end
      for i ∈ 1:σcnt
         qns = qnlab_full_simple(jmax,s,vtm)
         for j ∈ 1:size(vals,1)
            println(io, englin(s, vals[j,i]/csl, qns[j,:], i-1))
         end # j loop
      end # i loop 
   end # io
end

function linestrng(frql,qnu,σu,qnl,σl)   #this formats the lines for the transition writer
   part  = lpad(qnu[1]*0.5,4)*"," #J
   part *= lpad(qnu[2],3)*","     #N
   part *= lpad(qnu[3],3)*","     #Ka
   part *= lpad(qnu[4],3)*","     #Kc
   part *= lpad(qnu[5],3)*","     #vt
   part *= lpad(  σu-1,3)*","     #σ
   part *= lpad(qnl[1]*0.5,4)*"," #J
   part *= lpad(qnl[2],3)*","     #N
   part *= lpad(qnl[3],3)*","     #Ka
   part *= lpad(qnl[4],3)*","    #Kc
   part *= lpad(qnl[5],3)*","    #vt
   part *= lpad(  σl-1,3)*","     #σ
   part *= " "*@sprintf("%13.4f", frql[1])*","
   #part *= @sprintf("%10.4f", frql[4])*","
   part *= @sprintf("%12.6f", frql[2])*","
   part *= @sprintf("%10.4f", frql[3])
   return part
end

function writefreqs(molnam,ctrl,freqs,inds)
   qunus = qnlab_full_simple(ctrl.Jmax,ctrl.S,ctrl.vtmax)
   out = fill("0",size(freqs,1))
   for i in 1:size(freqs,1)               #writing the lines
      qunus[inds[i,1],:]
      inds[i,2]
      qunus[inds[i,3],:]
      inds[i,4]
      out[i] = linestrng(freqs[i,:], qunus[inds[i,1],:],inds[i,2], qunus[inds[i,3],:],inds[i,4])
   end
   io = open("$molnam.sim", "w") do io    #printing the lines
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.sim!")  #notice in terminal
end
