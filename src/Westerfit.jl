
module Westerfit

include("@__DIR__/../main.jl")

precompile(westerfit, String)

export westerfit, westersim, westereng, ctrlinp

"""
Hi! If you are trying to read the source code, I recommend you start in main.jl
It's broken up like this so for easier testing on my end as I just include
main.jl rather than calling the whole module when I need something quick.

The key files are:
userdefop.jl or "User Defineable Operators" has most of the Hamiltonian
common.jl is a collection of a lot of the smaller functions, especially the ones
   that relate to quantum numbers
assign.jl contains the routines for assigning quantum numbers after 
   diagonalization. It doesn't help with assigning QNs to the transitions
transitions.jl calculates the transition intensities & does the frequency filtering
optimizer.jl has the levenberg-marquadt implementation
files_in.jl and files_out.jl are for the input & output file handling respectively

Anything else is likely exerpimental or no longer used
"""

end # module
