
module Westerfit

include("@__DIR__/../main.jl")

precompile(westerfit, (String))
precompile(inp_reader, (String))
precompile(westereng, (String, Controls, Vector{Term}))
precompile(westersim, (String, Controls, Vector{MuOp}, Eigs)))
precompile(westerfit, (String, Controls, Vector{Term}))
#precompile(H_calc, (
#precompile(opt_calc, (

export westerfit, westersim, westereng, ctrlinp

"""
Hi! If you are trying to read the source code, I recommend you start in main.jl
It's broken up like this so for easier testing on my end as I just include
main.jl rather than calling the whole module when I need something quick.

The key files are:
main.jl has the big functions of westerfit (fitter), westereng (energy calculator)
   and westersim (simulator)
new_ham.jl has most of the Hamiltonian
common.jl is a collection of a lot of the smaller functions, especially the ones
   that relate to quantum numbers
assign.jl contains the routines for assigning quantum numbers after 
   diagonalization. It doesn't help with assigning QNs to the transitions
transitions.jl calculates the transition intensities & does the frequency filtering
optimizer.jl has the levenberg-marquadt implementation
files_in.jl and files_out.jl are for the input & output file handling respectively

Anything else is likely experimental or no longer used
"""

function julia_main()::Cint
  westerfit(ARGS)
  return 0 # if things finished successfully
end
Base.@ccallable function main()::Cint
   westerfit(ARGS)
   return 0
end

end # module
