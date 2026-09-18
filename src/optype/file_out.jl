
#### ENERGY FILE
function englin(s,eng,qunl,σ)
  # part = lpad(qunl[1],4)*"/2," #J
   part = lpad(qunl[1]*0.5,4)*","
   part *= lpad(qunl[2],4)*"," # N
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
   println(io, "   J,   N,  Ka,  Kc,  vt,   σ,"*lpad("Energy (cm-1)",17))
      if iszero(s)
      end
      for i ∈ 1:σcnt
         qns = qnlab_full_simple(jmax,s,vtm)
         for j ∈ 1:size(vals,1)
            println(io, englin(s, vals[j,i]/csl, qns[j,:], i))
         end # j loop
      end # i loop 
   end # io
end

#### SIMULATION FILE
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
   io = open(molnam*".sim", "w") do io    #printing the lines
   println(io, "  J', N',Ka',Kc',vt', σ',   J,  N, Ka, Kc, vt,  σ,"*
         lpad("ν (MHz),",15)*lpad("Int (arb),",13)*lpad("E_low (cm⁻¹)",14))
      for i in out
         println(io, i)
      end
   end
   println("Transitions written to $molnam.sim!")  #notice in terminal
end


##### RESIDUALS FILE
function reslin(qns,ofrq,omc,cfrq)       #this function does the pretty formatting for the
   part  = lpad(   (qns[ 1]),5)*"," # j     #reswriter function with a bunch of lpads
   part *= lpad(Int(qns[ 2]),3)*"," # n
   part *= lpad(Int(qns[ 3]),3)*"," # ka
   part *= lpad(Int(qns[ 4]),3)*"," # kc
   part *= lpad(Int(qns[ 5]),3)*"," # vt
   part *= lpad(Int(qns[ 6]),3)*"," # m
   part *= lpad(   (qns[ 7]),6)*","
   part *= lpad(Int(qns[ 8]),3)*","
   part *= lpad(Int(qns[ 9]),3)*","
   part *= lpad(Int(qns[10]),3)*","
   part *= lpad(Int(qns[11]),3)*","
   part *= lpad(Int(qns[12]),3)*","
   part *= " "*lpad(@sprintf("%0.5f", ofrq), 16)*"," #this part is where the number of sig figs gets fixed
   part *= " "*lpad(@sprintf("%0.6f", omc), 16)*","      #sometimes it breaks because of floating point errors
   part *= " "*lpad(@sprintf("%0.5f", cfrq), 16)
end
function reswritter(molnam::String,qns,lines::Lines,omcs,cfrqs)  #writes the .res file
   len = size(omcs,1)
   io = open("$molnam.res", "w") do io
      for i in 1:len
         println(io, reslin(qns[i,:], lines.frqs[i], omcs[i], cfrqs[i]))
      end                   
   end
end

#### OUTPUT FILE
function ctrl_print(io,ctrl::Controls)
   println(io, "Control Parameters")
   @inbounds for i ∈ propertynames(ctrl)
      println(io, i, " = ", getproperty(ctrl,i))
   end
end
function ham_init_print(io,ℋ::Vector{Term})
   perm = sortperm(map(x->x.nam, ℋ))
   for i ∈ perm
      println(io, lpad(ℋ[i].nam,12),";", lpad(ℋ[i].val, 30),";", lpad(ℋ[i].scl, 6))
   end
end
function sum_init_print(io, nparam, lins, omc)
   χ2 = sand(Diagonal(lins.wght)^2, omc)
   dof = length(omc) - nparam
   rms = √((omc' * omc)/dof)
   println(io, "\nNumber of lines = ",length(omc),", Number of parameters = ", nparam)
   println(io, "Initial χ2 = ", round(χ2, digits = 5),". The optimizer runs in terms of this value")
   println(io, "wrms = ", round(√(χ2/dof), digits=5), ", rms = ", round(rms,digits=5), " MHz\n")
end
function output_init(molnam::String,ctrl::Controls,ℋ::Vector{Term},lins)
   io = open(molnam*".out", "w") do io
   println(io, "westerfit!\n")
   println(io, molnam, "  @  ", now())
   println(io, "Sorry about the name...")
   ctrl_print(io, ctrl)
   println(io, "\nInitial Parameters")
   ham_init_print(io, ℋ)
   close(io)
   end
end
function ham_print(io,ℋ::Vector{Term},δ::Vector{Float64},prjct::Vector{Int})
   j = 1
   @inbounds for i ∈ eachindex(ℋ)
      if !iszero(ℋ[i].scl) && j ∉ prjct
         println(io, lpad(ℋ[i].nam,12),";", lpad(ℋ[i].val, 30),";", lpad(δ[j], 30))
         j += 1
      end
   end
end
function sum_print(io, iter, ℋ,δ, prjct, lnjct, rms, wrms, χ2)
   println(io, "\n",prod(fill("#",80)))
   println(io, "\nIteration: ", lpad(iter, 5))
   println(io, "χ2 = ", round(χ2, digits = 5))
   println(io, "wrms = ", round(wrms, digits=5), ", rms = ", round(rms,digits=5), " MHz\n")
   println(io, "Rejected Lines = ", length(lnjct)-1)
   if !iszero(length(lnjct)-1)
      println(io, lnjct)
   end
   println(io, "Temporarily frozen parameters = ", length(prjct)-1)
   if !iszero(prjct)
      println(io, getproperty.(ℋ[prjct] :nam))
   end
   println(io, "\nCurrent Parameter values:")
   ham_print(io,ℋ,δ, prjct)
end
function output_update(ℋ,δ,omc,w, counter)
   χ2 = sand(Diagonal(lins.wght)^2, omc)
   dof = length(omc) - nparam
   rms = √((omc' * omc)/dof)
   io = open(molnam*".out", "a") do io
      println(io, "After ",lpad(counter,3)," Iterations:\n")
      sum_print(io, ℋ, prjct, lnjct, rms, √(χ2/dof), χ2)
      println(io, "\n",repeat("-",30), "\n")
   close(io)
   end
end
function unc_arrange(ℋ::Vector{Term}, unc::Vector{Float64},prjct)#::Matrix{Float64}
   out = zeros(length(ℋ),3)
   j = 1
   deunit = unit_undict()
   for i ∈ 1:length(ℋ)
      if !iszero(ℋ[i].scl) j ∉ prjct
         out[i,1] = deunit[ℋ[i].unit](ℋ[i].val)
         out[i,2] = deunit[ℋ[i].unit](unc[j])
         out[i,3] = 1e3 * abs(unc[j] / ℋ[i].val)
         j += 1
      else
         out[i,1] = deunit[ℋ[i].unit](ℋ[i].val)
      end
   end
   perm = sortperm(getproperty.(ℋ, :nam))
   return perm, out[perm,:]
end
function ham_final_print(io, ℋ::Vector{Term}, unc, prjct)
   perm, out = unc_arrange(ℋ,unc,prjct)
   println(io, "Final parameter values & uncertanties & milicent\n")
   for i ∈ 1:length(perm)
      j = perm[i]
      if !iszero(ℋ[j].scl) j ∉ prjct
         mcnt =  @sprintf("%0.4f", out[i,3])
         println(io, lpad(ℋ[j].nam,12),"; ", lpad(out[i,1],30),"; ", lpad(out[i,2],30),"; ", lpad(mcnt,10))
      elseif !iszero(ℋ[j].scl) j ∉ prjct
         println(io, lpad(ℋ[j].nam,12),"; ", lpad(out[i,1],30),"; ", lpad("FROZEN BY CODE",30),
            "; ", lpad("UNDEFINED",10))
      else
         println(io, lpad(ℋ[j].nam,12),"; ", lpad(out[i,1],30),"; ", lpad("fixed",30),"; ", lpad("---",10))
      end
   end
end
function triangleprint(mat,nams;io=stdout,d=4,col=5)
   #io ≠ stdout ? io = open(io, "a") : io=stdout
   println(io,"\n\n  Correlation matrix:")
   l = size(mat,1)
   blocks = ceil(Int,l/col)
   for j in 1:blocks
   start = col*(j-1)+1
   stop = min(l,col*j)
   println(io,"\n"*prod(fill(" ",d))prod(lpad.(nams[start:stop],2d+2))*"\n")
   for i in start:l
      stop = min(i,col*j)
      part = lpad(nams[i],4)*prod(lpad.(round.(mat[start:stop,i],digits=d),2d+2))
      println(io,part)
   end; end
   io ≠ stdout ? println(io,"\n\n") : nothing
   #io ≠ stdout ? close(io) : nothing
end

function triangleprint(mat;io=stdout,d=4,col=5)
   nams = getproperty.(htrunc, :nam)
   triangleprint(mat,nams,io,d,col)
end

function output_final(molnam::String,ℋ,unc,corr, prjct, lnjct, rms, wrms, χ2, endp)
   io = open(molnam*".out", "a")
   println(io, "Fit completed by method of ", endp)
   # final rms
   println(io, "\n Final error values:")
   println(io, "χ2 = ", round(χ2, digits = 5))
   println(io, "wrms = ", round(rms, digits=5), ", rms = ", round(rms,digits=5), " MHz\n")
   # final pameters unc, millicent
   ham_final_print(io, ℋ, unc, prjct)
   # Journal formater
   # correlation matrix
   triangleprint(corr, getproperty.(ℋ[1:end .!= prjct], :nam), io=io)
   # apology
   println(io, "Again sorry about the name")
   close(io)
   println("output written to ", molnam,".out!")
end











