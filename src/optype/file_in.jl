
function inp_reader(molnam::String)
   inp = TOML.tryparsefile("$molnam.toml")
   if typeof(inp)==Base.TOML.ParserError
      @warn "Something is wrong with input TOML! Double check the file.
         Look at line $(inp.line) for a $(inp.type)"
      exit()
   end
   ctrl = controls_in(inp["controls"])
   #prms, scls, stgs 
#   secpart = secorder_in(inp["second_order"],length(ctrl.NFOLD))
   H = ops_in(inp["user_def"])
   μs = mus_in(inp["intensity"])
#   return ctrl, H, vcat(secpart[1], prms), vcat(secpart[2], scls), μs
   return ctrl, H, μs
end

function ctrl_sanity(ctrl::Controls)::Controls
   if ctrl.apology
      println("Sorry about the name...")
   end
   if length(ctrl.NFOLD) > 1
      @info "It appears you are using the n-top mode of westerfit.
      By doing so, you agree to not complain about the runtime to Julia (fka Wes).
      The matrices are very large & she's already losing sleep about it"
   end
   if iszero(ctrl.stages)
      println("You must define the number of stages")
      exit()
   end
   if iszero(length(ctrl.NFOLD))||(isone(length(ctrl.NFOLD))&&iszero(ctrl.NFOLD[1]))
      printstyled("Suppressing torsional behavior\n",color=:red)
      setproperty!(ctrl, :mcalc, 0)
      setproperty!(ctrl, :vtcalc, 0)
      setproperty!(ctrl, :vtmax, 0)
   end
   if isodd(2*ctrl.S)&&iseven(2*ctrl.Jmax)
      println("Jmax must be half integer for half interger S. Adding 1/2")
      setproperty!(ctrl, Jmax, ctrl.Jmax + 0.5)
   elseif iseven(2*ctrl.S)&&isodd(2*ctrl.Jmax)
      println("Jmax must be integer for interger S. Adding 1/2")
      setproperty!(ctrl, Jmax, ctrl.Jmax + 0.5)
   end
   if isone(ctrl.stages) && length(ctrl.NFOLD) > 1
      @warn "Single stage diagonalization is not implemented for multiple tops.
      The assignment process is non-sensical & the diagonalization times are unpleasant.
      One day I may try to implement this but I don't presently want to.
      This code will now terminate. Sorry"
      exit()
   end
   return ctrl
end
function controls_in(inp::Dict{String,Any})::Controls
   ctrl = Controls()
   for f ∈ eachindex(inp)
      setproperty!(ctrl, Symbol(f), inp[f])
   end
   return ctrl_sanity(ctrl)
end

#function secordinit_lim(topcount=0)::Dict{String,Int}
#   prd = Dict{String,Int}("A" => 1, "B" => 2, "C" => 3, "Dab" => 4,
#                          "Z" => 1, "X" => 2, "Y" => 3, "Dxz" => 4, 
#      "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵzx" => 8, "ϵxz" => 8,
#      "Czz" => 5, "Cxx" => 6, "Cyy" => 7, "Czx" => 8, "Cxz" => 8,
#      "χzz" => 9, "χxz" =>10, "χxmy"=>11, "χxx-χyy"=>11,
#        "α" => 9,   "δ" =>10,    "β"=>11)#,
#   if topcount ≥ 1
#      prd["F"]  = hccount + 1
#      prd["ρz"] = hccount + 2
#      prd["ρx"] = hccount + 3
#      prd["Vn"] = hccount + 4
#      prd["V3"] = hccount + 4
#   end
#   for i ∈ 1:topcount
#      prd["F_$i"]  = hccount + 1 + 4(i-1)
#      prd["ρz_$i"] = hccount + 2 + 4(i-1)
#      prd["ρx_$i"] = hccount + 3 + 4(i-1)
#      prd["Vn_$i"] = hccount + 4 + 4(i-1)
#      prd["V3_$i"] = hccount + 4 + 4(i-1)
#   end
#   return prd
#end
#function secnam_init()::Vector{String}
#   return ["AZ"; "BX"; "CY"; "Dab"; "ϵzz"; "ϵxx"; "ϵyy"; "ϵzx"; "χzz"; "χxx-χyy"; "χxz"]
#end
#function secorder_in(inp::Dict{String,Any},topcount=0)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int}}
#   params = zeros(hccount + 4*topcount)
#   scales = zeros(hccount + 4*topcount)
#   stages = zeros(Int, hccount + 4*topcount)
#   sodict = secordinit_lim(topcount)
#   for i ∈ eachindex(inp)
#      ind = sodict[i]
#      params[ind], scales[ind], stages[ind] = inp[i] 
#   end
#   return params, scales, stages
#end

function unit_dict()::Dict{String,Function}
   return Dict{String,Function}("MHz"=>x->x,
     "cm-1"=>x->csl*x,
      "kHz"=>x->1e-3*x,
       "Hz"=>x->1e-6*x,
      "mHz"=>x->1e-9*x,
      "GHz"=>x->1e3*x,
      "THz"=>x->1e6*x,
      "arb"=>x->x,
        "z"=>x->0.0,
       "eV"=>x->241_798_840.7662022*x, # <---- REPLACE WITH PHYS CONST
     "Hart"=>x->6_579_681_360.732768*x, # <---- REPLACE WITH PHYS CONST
       "nm"=>x->c_0.val * 1e3 / x,
       "μm"=>x->c_0.val/x,
         ""=>x->x)
end

function opfn_parse(x::String)
   q = [1;0]
   temp = split(x, r"\^|_")
   fn = getfield(Main, Symbol(temp[1]))
   if length(temp) > 3
      @warn "A user-defined operator got fucked up. Too many ^ & _
         Double check function input of $x"
   end
   for i ∈ 2:length(temp)
      q[i-1] = parse(Int,temp[i])
   end
   return fn,q
end
opfn_parse(x::SubString{String}) = opfn_parse(string(x))
#function opfn_proc(x,rf::Vector{OpFunc},tf::Vector{OpFunc})
function opfn_proc(x,rf,tf)#::Tuple(Vector{OpFunc},Vector{OpFunc})
   fn,q = opfn_parse(x)
   if methods(fn)[1].sig.parameters[2]==RPsi
      rf = vcat(rf, OpFunc(fn, q[1], q[2]) ) 
   elseif methods(fn)[1].sig.parameters[2]==TPsi
      if occursin("α",x)
         tf = vcat(tf, OpFunc(fn, q[1], 1) )
      elseif occursin("β",x)
         tf = vcat(tf, OpFunc(fn, q[1], 2) )
      elseif occursin("γ",x)
         tf = vcat(tf, OpFunc(fn, q[1], 3) )
      else
         tf = vcat(tf, OpFunc(fn, q[1], q[2]) )
      end
   else
      @warn "A function not dependent on the wavefunctions was called.
         Please double check the manual for valid functions.
         Name is $x"
   end
   return rf, tf
end
function op_parse(fstr::String,v::Float64)::Op
#   rf=Vector{OpFunc{Float64,RPsi}}[]
#   tf=Vector{OpFunc{Float64,TPsi}}[]
   rf = Union{}[]
   tf = Union{}[]
   xs = split(fstr, ' ')
   @inbounds for i ∈ eachindex(xs)
      rf,tf = opfn_proc(xs[i],rf,tf)
   end
   return Op(v,rf,tf)
end
#op_parse(fstr::SubString{String}) = op_parse(string(fstr))
function op_parse(fstr::Vector{String},vals::Vector{Float64})::Vector{Op}
   ops = Union{}[]
   @inbounds for i ∈ eachindex(fstr)
      ops = vcat(ops, op_parse(fstr[i],vals[i]))
   end
   return ops
end

prm_parse(v2::Float64)::Float64 = 1.0
function prm_parse(v2::Vector{Float64})::Vector{Float64}
   out = ones(length(v2))
   out[2:end] .= v2[2:end]
   return out
end

function ops_in(inp::Dict{String,Any})::Vector{Term}
   H = Vector{Term}(undef,length(inp))
   units = unit_dict()
   ln = 1
   for i ∈ eachindex(inp)
      vec = inp[i]
      if (typeof(vec[1]) == String)||typeof(vec[1]) == Vector{String}
         ops = op_parse(vec[1], prm_parse(vec[2]))
         H[ln] = Term(i, # name
            units[ vec[3] ](vec[2][1] ),
            ops, # operators
            vec[4], # scale
            vec[5]) # stage
      elseif typeof(vec[1]) == Vector{Any} && length(vec[1])==5
         opstr = vec[1][1]
         val = vec[1][2]
         for j ∈ 2:length(vec)
            opstr = vcat(opstr, vec[j][1])
            val = vcat(val, vec[j][2])
         end
         ops = op_parse(opstr, prm_parse(val))
         H[ln] = Term(i, # name
            units[ vec[1][3] ](vec[1][2] ),
            ops, # operator
            vec[1][4], # scale
            vec[1][5]) # stage
      else
         @warn "Something very weird happed with the input file at key $i"
      end # if
      ln += 1
   end
   return H
end

function mufn_proc(x,rf,tf)
   fn,q = opfn_parse(x)
   if methods(fn)[1].sig.parameters[2]==RPsi
      rf = vcat(rf, MuFunc(fn, q[1], q[2]) ) 
   elseif methods(fn)[1].sig.parameters[2]==TPsi
      if occursin("α",x)
         tf = vcat(tf, MuFunc(fn, q[1], 1) )
      elseif occursin("β",x)
         tf = vcat(tf, MuFunc(fn, q[1], 2) )
      elseif occursin("γ",x)
         tf = vcat(tf, MuFunc(fn, q[1], 3) )
      else
         tf = vcat(tf, MuFunc(fn, q[1], q[2]) )
      end
   else
      @warn "A function not dependent on the wavefunctions was called.
         Please double check the manual for valid functions.
         Name is $x"
   end
   return rf, tf
end
function mu_parse(fstr)
   rf=Vector{MuFunc}[]
   tf=Vector{MuFunc}[]
   xs = split(fstr, ' ')
   @inbounds for i ∈ eachindex(xs)
      rf,tf = mufn_proc(xs[i],rf,tf)
   end
   return rf,tf
end
function mus_in(inp::Dict{String,Any})
   μs = Vector{Mu}(undef,length(inp))
   ln = 1
   for i ∈ eachindex(inp)
      vec = inp[i]
      rf, tf = mu_parse(vec[1])
      μs[ln] = Mu(i, rf, tf, vec[2] )
      ln += 1
   end
   return μs
end
