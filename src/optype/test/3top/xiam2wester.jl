using DelimitedFiles

file = readdlm("26dma_x.lne")
out = zeros(size(file,1), 14)

function xs2w2(x::SubString)::Int
   if x=="S1"
      y = 0
   elseif x=="S2"
      y = 3
   elseif x=="S3"
      y = 2
   elseif x=="S4"
      y = 4
   elseif x=="S5"
      y = 5
   elseif x=="S6"
      y = 8
   elseif x=="S7"
      y = 12
   elseif x=="S8"
      y = 7
   elseif x=="S9"
      y = 11
   elseif x=="S10"
      y = 9
   else
      @show x
      y = -1
      println("FUCKFUCKFUCK")
   end
   return y
end

@show file[:,7]

###
out[:, 1] = file[:, 1] # J
out[:, 2] = file[:, 1] # N
out[:, 3] = file[:, 2] # Ka
out[:, 4] = file[:, 3] # Kc
out[:, 5] .= 0 # vt
out[:, 6] = xs2w2.(file[:, 7]) # σ
####
out[:, 7] = file[:, 4] # J
out[:, 8] = file[:, 4] # N
out[:, 9] = file[:, 5] # Ka
out[:,10] = file[:, 6] # Kc
out[:,11] .= 0 # vt
out[:,12] = xs2w2.(file[:, 7]) # σ
###
out[:,13] = file[:,end]
out[:,14] .= 0.002

writedlm("26dma.lne", out, ',')