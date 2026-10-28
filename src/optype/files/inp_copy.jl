function iterativecopier(molnam::String)   #this iterates the files up a value
   oldarray = ["8", "7", "6", "5", "4", "3", "2", "1", ""]
   newarray = ["9", "8", "7", "6", "5", "4", "3", "2", "1"]
   for i in 1:length(oldarray)
      filename = string(molnam, oldarray[i], ".toml")
      newfilename = string(molnam, newarray[i], ".toml")
      if isfile(filename) == true
         cp(filename, newfilename, force=true)
      end
   end
end


function inp_info_print(io,molnam,inf,rms,wrms,numlin)
   println(io,"[info]")
   inf["molnam"] = molnam
   inf["rms"] = rms
   inf["wrms"] = wrms
   inf["line_count"] = numlin
   TOML.print(io, inf)
end

function ctrl2strng(ctrl,i)
   x = getproperty(ctrl,i)
   if typeof(x)==String
      x = "\""*x*"\""
   end
   return x
end
function inp_ctrl_print(io,ctrl::Controls)
   println(io, "\n[controls]")
   for i ∈ propertynames(ctrl)
      println(io, "\"", i, "\" = ", ctrl2strng(ctrl,i) )
   end
end
function opf2strng(o)::String
   out = string(o.f.obj.x)
   if !iszero(o.q)
      out *= "^"*string(o.l)*"_"*string(o.q)
   elseif o.l > 1 && iszero(o.q)
      out *= "^"*string(o.l)
   else
   end
   return out 
end
function op2strng(O)::String
   out = "\""
   temp = vcat(O.rf,O.tf)
   for i ∈ 1:length(temp)-1
      out *= opf2strng(temp[i])*" "
   end
   out *= opf2strng(temp[end])
   return out*"\""
end
function op2strng(O::Vector)::String
   if isone(length(O))
      out = ""
      out *= op2strng(O[1])
   else
      out = "["
      for i ∈ 1:length(O)-1
         out *= op2strng(O[i])*", "
      end
      out *= op2strng(O[end])*"]"
   end
   return out
end
function val2strng(T::Term,udict)::String
   if isone(length(T.ops))
      out = string(udict[T.unit](T.val))
   else
      out = "["*string(udict[T.unit](T.val))
      for i ∈ 2:length(T.ops)
         out*= ", "*string(T.ops[i].val)
      end
      out *= "]"
   end
   return out
end

function inp_ham_print(io,ℋ)
   udict = unit_undict()
   println(io, "\n[hamiltonian]")
   println(io, "#Nam = [Func, val, unit, scale, stage]")
   println(io, "# α β γ δ Δ ϵ η λ Λ μ ϕ Φ σ ω")
   for i ∈ sortperm(getproperty.(ℋ, :nam))
      println(io, "\"", ℋ[i].nam, "\" = [",op2strng(ℋ[i].ops),", ",
         val2strng(ℋ[i], udict),", \"", ℋ[i].unit,"\", ", ℋ[i].scl,", ", ℋ[i].stg,"]")
   end
end

function inp_mu_print(io,μs)
   if μs==nothing
      println(io, "\n[intensity]")
   else
      println(io, "\n[intensity]")
      for i ∈ sortperm(getproperty.(μs, :nam))
         println(io, "\"", μs[i].nam, "\" = [", op2strng(μs[i]),", ", μs[i].val,"]")
      end
   end
end

function inp_print(io,molnam::String,info::Dict,rms::Float64,wrms::Float64,numlin::Int,
                     ctrl::Controls,ℋ::Vector{Term},μs)
   inp_info_print(io,molnam,info,rms,wrms,numlin)
   inp_ctrl_print(io,ctrl)
   inp_ham_print(io,ℋ)
   inp_mu_print(io,μs)
end