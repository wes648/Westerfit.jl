@kwdef mutable struct Controls
   apology::Bool = true
   RUNmode::String = "ESF"
   stages::Int = 0
   Irrep::String = "Ir"
   assign::String = "ram36"
   NFOLD::Vector{Int} = [0] # vector of symmetry folds of rotors
   S::Float64 = 0. # float for spin value could maybe turn into int for 2s but eh
   Jmax::Float64 = 0. # maximum J value
   mcalc::Int = 10 # maximum |m| value for free rotor basis, basis size will be 2mmax+1
   vtcalc::Int = 8 # maximum vt state output by second diag stage & to be used in third. basis size will be vtmax+1
   vtmax::Int = 0 # maximum vt state output by final diagonalization stage. basis size will be vtmax+1
   νmin::Float64 = 0.2
   νmax::Float64 = 40.0
   TK::Float64 = 8. # temperature in Kelvin to be used in simulation
   INTthres::Float64 = 0.0001
   λlm0::Float64 = 0.001
   turducken::Int = 1
   maxiter::Int = 60
   BOLD::Int = 0
   REJECT::Float64 = 10.0
   goal::Float64 = 1.0
   overwrite::Bool = true 
end

function inp_reader(molnam::String)
   inp = TOML.tryparsefile("$molnam.toml")
   if typeof(inp)==Base.TOML.ParserError
      @warn "Something is wrong with input TOML! Double check the file.
         Look at line $(inp.line) for a $(inp.type)"
      exit()
   end
   ctrl = controls_in(inp["controls"])
   #prms, scls, stgs 
   secpart = secorder_in(inp["second_order"])
   H, prms, scls = ops_in(inp["user_def"])
   return ctrl, H, vcat(secpart[1], prms), vcat(secpart[2], scls)
end

function ctrl_sanity(ctrl::Controls)::Controls
   if iszero(ctrl.stages)
      println("You must define the number of stages")
      exit()
   end
   if iszero(length(ctrl.NFOLD))||(isone(length(ctrl.NFOLD))&&iszero(ctrl.NFOLD[1]))
      println("Suppressing torsional behavior")
      setproperty!(ctrl, mcalc, 0)
      setproperty!(ctrl, vtcalc, 0)
      setproperty!(ctrl, vtmax, 0)
   end
   if ctrl.νmin > ctrl.νmax
      println("Enforcing νmin < νmax")
      temp = ctrl.νmin
      setproperty!(ctrl, νmin, ctrl.νmax)
      setproperty!(ctrl, νmax, temp)
   end
   if isodd(2*ctrl.S)&&iseven(2*ctrl.Jmax)
      println("Jmax must be half integer for half interger S. Adding 1/2")
      setproperty!(ctrl, Jmax, ctrl.Jmax + 0.5)
   elseif iseven(2*ctrl.S)&&isodd(2*ctrl.Jmax)
      println("Jmax must be integer for interger S. Adding 1/2")
      setproperty!(ctrl, Jmax, ctrl.Jmax + 0.5)
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

function secordinit_lim(topcount=0)::Dict{String,Int}
   prd = Dict{String,Int}("A" => 1, "B" => 2, "C" => 3, "Dab" => 4,
                          "Z" => 1, "X" => 2, "Y" => 3, "Dxz" => 4, 
      "ϵzz" => 5, "ϵxx" => 6, "ϵyy" => 7, "ϵzx" => 8, "ϵxz" => 8,
      "Czz" => 5, "Cxx" => 6, "Cyy" => 7, "Czx" => 8, "Cxz" => 8,
      "χzz" => 9, "χxz" =>10, "χxmy"=>11, "χxx-χyy"=>11,
        "α" => 9,   "δ" =>10,    "β"=>11)#,
   if topcount ≥ 1
      prd["F"]  = hccount + 1
      prd["ρz"] = hccount + 2
      prd["ρx"] = hccount + 3
      prd["Vn"] = hccount + 4
   end
   for i ∈ 1:topcount
      prd["F_$i"]  = hccount + 1 + 4(i-1)
      prd["ρz_$i"] = hccount + 2 + 4(i-1)
      prd["ρx_$i"] = hccount + 3 + 4(i-1)
      prd["Vn_$i"] = hccount + 4 + 4(i-1)
   end
   return prd
end
function secnam_init()::Vector{String}
   return ["AZ"; "BX"; "CY"; "Dab"; "ϵzz"; "ϵxx"; "ϵyy"; "ϵzx"; "χzz"; "χxx-χyy"; "χxz"]
end
function secorder_in(inp::Dict{String,Any},topcount=0)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Int}}
   params = zeros(hccount + 4*topcount)
   scales = zeros(hccount + 4*topcount)
   stages = zeros(Int, hccount + 4*topcount)
   sodict = secordinit_lim(topcount)
   for i ∈ eachindex(inp)
      ind = sodict[i]
      params[ind], scales[ind], stages[ind] = inp[i] 
   end
   return params, scales, stages
end

function unit_dict()::Dict{String,Float64}
   return Dict{String,Float64}("MHz"=>1.,"cm-1"=>29979.2458,"kHz"=>1e-3,"Hz"=>1e-6,
   "mHz"=>1e-9,"GHz"=>1e3,"THz"=>1e6,"arb"=>1.,"z"=>0.0,
   "eV"=>241_798_840.7662022,"Hart"=>6_579_681_360.732768)
end

function opfn_parse(x)
   q = [1;0]
   temp = split(x, r"\^|_")
   fn = getfield(Main, Symbol(temp[1]))
   if length(temp) > 3
      @warn "A user-defined operator got fucked up. Too many ^ & _
         Double check function input"
   end
   for i ∈ 2:length(temp)
      q[i-1] = parse(Int,temp[i])
   end
   return fn,q
end
function opfn_proc(x,rf,tf)
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
         Please double check the manual for valid functions"
   end
   return rf, tf
end
function op_parse(fstr)
   rf=Vector{OpFunc}[]
   tf=Vector{OpFunc}[]
   xs = split(fstr, ' ')
   @inbounds for i ∈ eachindex(xs)
      rf,tf = opfn_proc(xs[i],rf,tf)
   end
   return rf,tf
end

function ops_in(inp::Dict{String,Any})
   H = Vector{Op}(undef,length(inp))
   prm = zeros(Float64, length(inp))
   scl = zeros(Float64, length(inp))
   units = unit_dict()
   ln = 1
   for i ∈ eachindex(inp)
      vec = inp[i]
      prm[ln] = vec[2] * units[ vec[3] ]
      rf,tf = op_parse(vec[1])
      H[ln] = Op(i, rf, tf, vec[5])
      scl[ln] = vec[4]
      ln += 1
   end
   return H, prm, scl
end
