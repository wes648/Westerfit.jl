
"""
This file is the actual implementation of the information in types.jl & baseops.jl
"""

using DelimitedFiles

include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../file_in.jl")

const csl = 29979.2458

function westereng(molnam::String)
   ctrl = ctrlinp(molnam)
   val, errs,ℋ, stgs = secordinp(molnam,ctrl)
   ℋ,stgs, errs = opreader(molnam,ctrl,ℋ,stgs,errs)
   ℋ = stgvalset(ℋ,stgs)
end