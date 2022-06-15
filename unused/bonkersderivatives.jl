function dHdA(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   B = rp[2]
   sδ = sin(rp[4])
   s2δ = sin(2.0*rp[4])
   cδ = cos(rp[4])
   F = rp[5]
   #dAeff/dAp
   drp[1] = ((B^2)*(F^2)*(cδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dBr/dAp
   drp[2] = B / (A*(cδ^2) + B*(sδ^2)) - ((A*B) / ((A*(cδ^2) + B*(sδ^2))^2))*(cδ^2)
   #dDab/dAp
   drp[4] =((A^2)*(cδ^2)*(sδ^2)*s2δ + 2.0A*B*(cδ^4)*s2δ + 2.0A*B*(sδ^4)*s2δ +
      3.0(B^2)*(cδ^2)*(sδ^2)*s2δ - (A^2)*(cδ^4)*s2δ - (A^2)*(sδ^4)*s2δ -
      2.0A*B*(cδ^2)*(sδ^2)*s2δ) / (4(A^4)*(B^2))
   #dFr/dAp
   drp[5] = ((F^2)*(B*(cδ^2) + A*(sδ^2))) / (B*F*(cδ^2) + A*F*(sδ^2) - A*B)
   #dFρ/dAp
   drp[6] = ((B^2)*(F^2)*(cδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out[1]
end
function dHdB(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   B = rp[2]
   sδ = sin(rp[4])
   s2δ = sin(2.0*rp[4])
   cδ = cos(rp[4])
   F = rp[5]
   #dAeff/dBp
   drp[1] = ((A^2)*(F^2)*(sδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dBr/dBp
   drp[2] = ((A^2)*(cδ^2)) / ((A*(cδ^2) + B*(sδ^2))^2)
   #dDab/dBp
   drp[4] =((B^2)*(cδ^4)*s2δ + (B^2)*(sδ^4)*s2δ + 2.0A*B*(cδ^2)*(sδ^2)*s2δ -
      (B^2)*(cδ^2)*(sδ^2)*s2δ - 2.0A*B*(cδ^4)*s2δ - 2.0A*B*(sδ^4)*s2δ -
      3.0(A^2)*(cδ^2)*(sδ^2)*s2δ) / (4(A^2)*(B^4))
   #dFr/dBp
   drp[5] = ((A^2)*(F^2)*(sδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dFρ/dBp
   drp[6] = ((A^2)*(F^2)*(sδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out[1]
end
function dHdδ(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   B = rp[2]
   sδ = sin(rp[4])
   s2δ = sin(2.0*rp[4])
   cδ = cos(rp[4])
   c2δ = cos(2.0*rp[4])
   F = rp[5]
   #dAeff/dδ
   drp[1] = (A*B*F*(B*F*s2δ - A*F*s2δ)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dBr/dδ
   drp[2] = (A*B*(A*s2δ - B*s2δ)) / ((A*(cδ^2) + B*(sδ^2))^2)
   #dDab/dδ
   drp[4] = ((A^3)*(cδ^2)*(s2δ^2) + (B^3)*(sδ^2)*(s2δ^2) + 2.0*B*(A^2)*(cδ^4)*c2δ
      + 2.0*B*(A^2)*(sδ^4)*c2δ + 3.0*A*(B^2)*(cδ^2)*(s2δ^2) + 3.0*B*(A^2)*(sδ^2)*(s2δ^2)
      + 2.0*(A^3)*(cδ^2)*(sδ^2)*c2δ + 2.0*A*(B^2)*(cδ^2)*(sδ^2)*c2δ
      - (B^3)*(cδ^2)*(s2δ^2) - (A^3)*(sδ^2)*(s2δ^2) - 3.0*A*(B^2)*(sδ^2)*(s2δ^2)
      - 3.0*B*(A^2)*(cδ^2)*(s2δ^2) - 2.0*A*(B^2)*(cδ^4)*c2δ - 2.0*(B^3)*(cδ^2)*(sδ^2)*c2δ
      - 2.0*A*(B^2)*(sδ^4)*c2δ - 2.0*B*(A^2)*(cδ^2)*(sδ^2)*c2δ) / (4(A^3)*(B^3))
   #dFr/dδ
   drp[5] = (A*(B^2)*(F^2)*s2δ - B*(A^2)*(F^2)*s2δ) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dFρ/dδ
   drp[6] = (A*B*F*(B*F*s2δ - A*F*s2δ)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out[1]
end
function dHdF(j,s,σ,vec,rp)
   drp = zeros(Float64, length(rp)+1)
   A = rp[1]
   B = rp[2]
   sδ = sin(rp[4])
   s2δ = sin(2.0*rp[4])
   cδ = cos(rp[4])
   c2δ = cos(2.0*rp[4])
   F = rp[5]
   #dAeff/dF
   drp[1] = (-(A^2)*(B^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dFr/dF
   drp[5] = ((B^2)*(F^2)*(cδ^4) + (A^2)*(F^2)*(sδ^4) + 2.0A*B*(F^2)*(cδ^2)*(sδ^2)
   - 2.0A*F*(B^2)*(cδ^2) - 2.0B*F*(A^2)*(sδ^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   #dFρ/dF
   drp[6] = (-(A^2)*(B^2)) / ((B*F*(cδ^2) + A*F*(sδ^2) - A*B)^2)
   U = ur(j,s,mcalc)
   if σ==0
      U *= ut(mcalc,j,s)
   end
   mat = Matrix(U*Htsr(drp,j,s,mcalc,σ)*U)
   out = transpose(vec)*mat*vec
   return out[1]
end
#interal parameters:
#      1  2  3   4    5   6   7   8  9 10 11 12
#prs =[A; B; C; Dab; Fr; Fρ; V3; ao; a; b; d; η]
#input parameters:
#      1  2  3  4  5   6   7    8    9    10 11
#prs =[A; B; C; δ; F; V3; ϵzz; ϵxx; ϵyy; ϵxz; η]
function build_jcbn_pam(inds,vecs,params)
   #fitprm = params[perm]
   jcbn = zeros(Float64,size(inds)[1],length(params))
   for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      #println(size(vecu))
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      #println(vecs[:,inds[a,3],σu+1])
      #println(transpose(vecl)*vecl)
      jcbn[a,1] = dHdA(ju,S,σu,vecu,params) - dHdA(jl,S,σl,vecl,params)
      jcbn[a,2] = dHdB(ju,S,σu,vecu,params) - dHdB(jl,S,σl,vecl,params)
      jcbn[a,3] = anaderiv(ju,S,σu,vecu,params,3) - anaderiv(jl,S,σl,vecl,params,3)
      jcbn[a,4] = dHdδ(ju,S,σu,vecu,params) - dHdδ(jl,S,σl,vecl,params)
      jcbn[a,5] = dHdF(ju,S,σu,vecu,params) - dHdF(jl,S,σl,vecl,params)
      for b in 6:length(params)
         jcbn[a,b] = anaderiv(ju,S,σu,vecu,params,b+1) - anaderiv(jl,S,σl,vecl,params,b+1)
      end
   end
   return jcbn
end
