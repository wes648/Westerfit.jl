
using DelimitedFiles

function pred2lne(sim::Array{Float64,2},s::Number)
"""
Converts the transitions output from westersim into the line format for
   westerfit. Allows for quick test scripting
"""
   if s != zero(s)
      out = zeros(Float64,size(sim)[1],12)
      out[:,1:10] = sim[:,1:10]
      out[:,11] = sim[:,11]
      out[:,12] = fill(0.08,size(sim)[1])
   else#this should be though
      out = zeros(Float64,size(sim)[1],12)
      out[:,1] = sim[:,1]
      out[:,2:5] = sim[:,1:4]
      out[:,6] = sim[:,5]
      out[:,7:10] = sim[:,5:8]
      out[:,11] = sim[:,9]
      out[:,12] = fill(0.08,size(sim)[1])
   end
   return out
end

function pred2lne(nam::String)
   f = readdlm("$nam.sim",',')
   writedlm("test_input.lne", pred2lne(f,0.0), ',')
end

pred2lne("2ba")
