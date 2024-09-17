#function quickfreqs(vals,inds,ofreqs)
#   cfreqs = zero(ofreqs)
#   @threads for i in 1:size(cfreqs,1)
#      cfreqs[i] = vals[inds[i,3],inds[i,2]+1] - vals[inds[i,6],inds[i,5]+1]
#   end
#   return cfreqs
#end
#construct jacobian
#function anaderiv(j,s,σ,vec,rp,rpid)
#   rp = zero(rp)
#   rp[rpid] = 1.0
#   U = ur(j,s,mcalc)
#   if σ==zero(σ)
#      U *= ut(mcalc,j,s)
#   end
#   mat = Matrix(U*Htsrmat(rp,j,s,mcalc,σ)*U)
#   out = transpose(vec)*mat*vec
#   return out
#end
