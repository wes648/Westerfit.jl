
function englin(s,eng,qunl)
   if s==zero(s)
      part = lpad(qunl[2],4)*"," # N
   else
     # part = lpad(qunl[1],4)*"/2,"
      part = lpad(qunl[1],4)*","  # 2J
      part *= lpad(qunl[2],4)*"," # N
   end
   part *= lpad(qunl[3],4)*"," # Ka
   part *= lpad(qunl[4],4)*"," # Kc
   part *= lpad(qunl[5],4)*"," #  m
   part *= lpad(qunl[6],4)*"," #  σ
####0.10f is what BELGI uses, 0.6f is for spcat
   part *= " "*lpad(@sprintf("%0.10f", eng), 16) # energy
   return part
end

"""
Outputs energy levels with state assignments to a csv-like file
"""
function engwriter(molnam,ctrl,energies,qunus)
   eng = energies[:,1]
   qns = qunus[:,:,1]
   for sc in 2:σcount(ctrl["NFOLD"])
      eng = vcat(eng,energies[:,sc])
      qns = vcat(qns,qunus[:,:,sc])
   end
   len = size(eng,1)
   out = fill("0",len)
   for i in 1:len
      energy = eng[i]/csl
      out[i] = englin(ctrl["S"],energy,qns[i,:])
   end
   io = open("$molnam.eng", "w") do io
      for i in out
         println(io, i)
      end
   end
end
