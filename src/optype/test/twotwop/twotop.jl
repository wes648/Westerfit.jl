
using LinearAlgebra

function htor(f,v,mc,σ)
   ms = (-3*mc+σ):3:(3*mc+σ)
   h = diagm(0=> f .*(ms .^2) .+ 0.5*v, 
      -1=>fill(-0.25v,length(ms)-1), 1=>fill(-0.25v,length(ms)-1))
   return h
end

function sin(mc)
   l = 2mc
   out = diagm(1=>fill(0.5im, l),-1=>fill(-0.5im, l))
   return out
end

function twotop(σ1,σ2)
   mc = 8
   f1 = 5.554669
   f2 = 5.523464
   v31 = 101.0
   v32 = 422.0
   f12 = 0.664
   vcc = -7.0
   vss = 68.0
   l = 2mc +1
   out = kron(I(l), htor(f1,v31,mc,σ1)) + kron(htor(f2,v32,mc,σ2), I(l))
   out += Diagonal(kron(f12, (-3*mc+σ2):3:(3*mc+σ2), (-3*mc+σ1):3:(3*mc+σ1)))
   out += kron(htor(0,vcc,mc,σ1), htor(0,1,mc,σ2))
   out += kron(vss, sin(mc), sin(mc))
   out *= 0.5
   out += out'
   #return eigvals(out)[1]
   return out
end

function rot()
   A = 0.3818404
   B = 0.1401459
   C = 0.1025080
   out = [A+0.5(B+C) 0 0.5(B-C); 0 B+C 0; 0.5*(B-C) 0 A+0.5(B+C)]
   return out
end

function hrt(σ1,σ2)
   mc = 8
   Q1 = -0.68129
   Q2 = -0.70434
   e0 = eigvals(twotop(σ1,σ2))[1]
   h = kron(twotop(σ1,σ2), I(3)) + kron(I((2mc+1)^2), rot())
   h += Diagonal(kron(Q1, ones(2mc+1), (-3*mc+σ1):3:(3*mc+σ1), -1:1))
   h += Diagonal(kron(Q2, (-3*mc+σ2):3:(3*mc+σ2), ones(2mc+1), -1:1))
   h *= 0.5
   h += h'
   return [e0; eigvals(h)[1:3]]
end