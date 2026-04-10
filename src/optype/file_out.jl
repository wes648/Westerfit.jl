

function englin(s,eng,qunl)
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
   part *= lpad(qunl[6],4)*"," # σ
   part *= " "*lpad(@sprintf("%0.10f", eng), 16) # eng
   return part
end

function engwriter(molnam, jmax,s,vtm,vals)
   σcnt = size(vals,2)
   io = open("$molnam.eng", "w") do io
      if iszero(s)
      end
      for i ∈ 1:σcnt
         qns = qnlab_full_simple(jmax,s,vtm,i-1)
         for j ∈ 1:size(vals,1)
            println(io, englin(s, vals[j,i]/csl, qns[j,:]))   #prints the line in the energy file
         end # j loop
      end # i loop 
   end # io
end