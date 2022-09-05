using LinearAlgebra, SparseArrays, Symbolics
using Base.Threads

#borrowed from common.jl to have an isolated test file
function Δlist(J,S)
"""
This generates a list of all J₃ values that satisfy the Triangle Conditions with
   inputs J₁ and J₂. This is typically used to determine N values for a given
   J and S.
"""
   max = Int64(J+S)
   min = Int64(abs(J-S))
   return collect(min:max)
end
function srprep(J,S)
"""
This calculates a lot of values helpful for the spin-rotation matrix structure:
   the list of n values, the number of K states for each n, an array of the
   start and end indices for each n in the list, and the number of substates for
   the given J and S pair. The array of indices has starts on the [:,1] values
   and the ends on the [:,2] values.
"""
   ns = Δlist(J,S)
   nd = 2 .* Int.(ns) .+ 1
   out = ones(Int, length(ns),2)
   out[1,2] = nd[1]
   for i in 2:length(ns)
      out[i,1] = out[i-1,2] + 1
      out[i,2] = out[i,1] + nd[i] - 1
   end
   jd = Int((2*S+1)*(2*J+1))
   return ns, nd, out, jd
end

#associated functions
function C0(F,I,J)
   out = @. F*(F+1) - I*(I+1) - J*(J+1)
end
function C1(F,I,J)
   out = C0(F,I,J)
   out = @. 3*out*(out+1.0)/4.0 - I*(I+1.0)*J*(J+1.0)
   out = @. out/(2.0*I*(2.0*I - 1.0)*J*(J+1.0)*(2.0*J-1.0)*(2.0*J+3.0))
end
function C2(F,I,J)
   if J==zero(J)
      return 0.0
   else
   out = (F*(F+1.0)-I*(I+1.0) - J*(J+2.0))/(8.0*I*(2.0*I-1.0)*J*(J+1.0)*(J+2.0))
   out = out*√(((F+I+J+2.0)*(J+I-F+1.0)*(F+J-I+1.0)*(F+I-J))/((2.0*J+1.0)*(2.0*J+3.0)))
   return out
end
end
function C3(F,I,J)
   out = √(((F+I-J)*(F+I+J+2.0)*(F+I+J+3.0))/((2.0*J+1.0)*(2.0*J+5.0)))
   out = out*(16.0*I*(2.0*I-1.0)*(J+1.0)*(J+2.0)*(2.0*J+3.0))^-1
   out = out*√((I-F+J+1.0)*(I-F+J+2.0)*(F-I+J+1.0)*(F-I+J+2.0)*(F+I-J-1.0))
   return out
end

#matrix elements
function Hhyp0K0J(pr,c1,J,K)
# <J K I F| H_hyp |J K I F>
   out = @. c1 * (3.0*K^2 - J*(J+1.0))*pr[13]
end
function Hhyp1K0J(pr,c1,J,K)
# <J K+1 I F| H_hyp |J K I F>
   out = @. -c1*(2.0*K+1.0)*√(1.5*(J+K+1.0)*(J-K))*pr[15]
# chi_=-1
end
function Hhyp2K0J(pr,c1,J,K)
# <J K+2 I F| H_hyp |J K I F>
   out = @. c1*√(1.5*(J+K+1.0)*(J+K+2.0)*(J-K)*(J-K-1.0))*pr[14]
end
function Hhyp0K1J(pr,c2,J,K)
# <J+1 K I F| H_hyp |J K I F>
   out = @. -3.0*c2*K*√((J+1.0)^2 - K^2)*pr[13]
end
function Hhyp1K1J(pr,c2,J,K)
# <J+1 K+1 I F| H_hyp |J K I F>
   out = @. -c2*(J - 2.0*K)*√(1.5*(J+K+2.0)*(J+K+1.0))*pr[15]
end
function Hhypm1K1J(pr,c2,J,K)
# <J+1 K-1 I F| H_hyp |J K I F>
   out = @. -c2*(J + 2.0*K)*√(1.5*(J-K+2.0)*(J-K+1.0))*pr[15]
end
function Hhyp2K1J(pr,c2,J,K)
# <J+1 K+2 I F| H_hyp| J K I F>
   out = @. c2*√(1.5*(J-K)*(J+K+1.0)*(J+K+2.0)*(J+K+3.0))*pr[14]
end
function Hhypm2K1J(pr,c2,J,K)
# <J+1 K-2 I F| H_hyp |J K I F>
   out = @. -c2*√(1.5*(J+K)*(J-K+1.0)*(J-K+2.0)*(J-K+3.0))*pr[14]
end
function Hhyp0K2J(pr,c3,J,K)
# <J+2 K I F| H_hyp |J K I F>
   out = @. 3*c3*√(((J+1.0)^2-K^2)*((J+2.0)^2 - K^2))*pr[13]
end
function Hhyp1K2J(pr,c3,J,K)
# <J+2 K+1 I F| H_hyp |J K I F>
   out = @. -2.0*c3*√(1.5*(J-K+1.0)*(J+K+1.0)*(J+K+2.0)*(J+K+3.0))*pr[15]
end
function Hhypm1K2J(pr,c3,J,K)
# <J+2 K-1 I F| H_hyp |J K I F>
   out = @. 2.0*c3*√(1.5*(J-K+1.0)*(J+K+1.0)*(J-K+2.0)*(J-K+3.0))*pr[15]
end
function Hhyp2K2J(pr,c3,J,K)
# <J+2 K+2 I F| H_hyp |J K I F>
   out = @. c3*√(1.5*(J+K+1.0)*(J+K+2.0)*(J+K+3.0)*(J+K+4.0))*pr[14]
end
function Hhypm2K2J(pr,c3,J,K)
# <J+2 K-2 I F| H_hyp |J K I F>
   out = @. c3*√(1.5*(J-K+1.0)*(J-K+2.0)*(J-K+3.0)*(J-K+4.0))*pr[14]
end


#matrix builder

function Hhyp0J(pr,F,I,J)
   c1 = C1(F,I,J)
   if J==zero(J)
      Hypmat = spzeros(typeof(pr[13]),1,1)
   else
      karray = collect(Float64,-J:J)
      ondiags = Hhyp0K0J(pr,c1,J,karray)
      of1diag = Hhyp1K0J(pr,c1,J,karray[1:end-1])
      of2diag = Hhyp2K0J(pr,c1,J,karray[1:end-2])
      Hypmat = spdiagm(length(karray),length(karray),0=>ondiags,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
   end
   return Hypmat
end
function Hhyp1J(pr,F,I,J)
   if J==zero(J)
      return spzeros(1,3)
   else
   c2 = C2(F,I,J)
   karray = collect(Float64,-J:J)
   pm2 = Hhypm2K1J(pr,c2,J,karray[2:end])
   pm1 = Hhypm1K1J(pr,c2,J,karray)
   p0 = Hhyp0K1J(pr,c2,J,karray)
   p1 = Hhyp1K1J(pr,c2,J,karray)
   p2 = Hhyp2K1J(pr,c2,J,karray[1:end-1])
   Hypmat = spdiagm(2*Int(J)+1,2*Int(J)+3, -1=>pm2, 0=>pm1, 1=>p0, 2=>p1, 3=>p2)
   return Hypmat
   end
end
function Hhyp2J(pr,F,I,J)
   c3 = C3(F,I,J)
   karray = collect(Float64,-J:J)
   pm2 = Hhypm2K2J(pr,c3,J,karray)
   pm1 = Hhypm1K2J(pr,c3,J,karray)
   p0 = Hhyp0K2J(pr,c3,J,karray)
   p1 = Hhyp1K2J(pr,c3,J,karray)
   p2 = Hhyp2K2J(pr,c3,J,karray)
   #I don't know if thsi next line is correct
   Hypmat = spdiagm(2*Int(J)+1,2*Int(J)+5, 0=>pm2, 1=>pm1, 2=>p0, 3=>p1, 4=>p2)
   return Hypmat
end

function Hhyp!(pr,mat,J,S)
   if S ≥ 1
   ns, nd, ni, jd = srprep(J,S)
   mat[1:nd[1],1:nd[1]] .+= Hhyp0J(pr,J,S,ns[1])
   for i in 2:length(ns)
      n = ns[i]
      n1part = Hhyp1J(pr,J,S,n-1.0)
      mat[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]]  .+= n1part
      mat[    ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]]  .+= Hhyp0J(pr,J,S,n)
      mat[    ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] .+= transpose(n1part)
      if i ≥ 3
         println(i)
         n2part = Hhyp2J(pr,J,S,n-2.0)
         mat[ni[i-2,1]:ni[i-2,2],   ni[i,1]:ni[i,2]] .+= n2part
         mat[   ni[i,1]:ni[i,2],ni[i-2,1]:ni[i-2,2]] .+= transpose(n2part)
      end
   end
   end
   return mat
end

#test section
@variables χ0 χ1 χ2
params = zeros(Number,15)
params[13:15] = [χ0, χ2, χ1]

out = spzeros(typeof(params[13]),15,15)
out = Hhyp!(params,out,2,1)
