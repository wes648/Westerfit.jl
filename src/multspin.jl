using LinearAlgebra, SparseArrays, WIGXJPFjl, BenchmarkTools
jnred(j::Real,n::Real)::Float64 = √((2*j+1)*(2*n+1))
nred(n::Real)::Float64 = √(n*(n+1)*(2*n+1))
function Δlist(J,S)
   return collect(Int(abs(J-S)):Int(J+S))
end
function srprep(J,S)
   ns = Δlist(J,S)
   nd = 2 .* Int.(ns) .+ 1
   ni = ones(Int, length(ns),2)
   ni[1,2] = nd[1]
   for i in 2:length(ns)
      ni[i,1] = ni[i-1,2] + 1
      ni[i,2] = ni[i,1] + nd[i] - 1
   end
   jd = Int((2.0*S+1.0)*(2.0*J+1.0))
   return ns, nd, ni, jd
end
function qngen(j,s)
   ns, nd, ni, jsd = srprep(j,s)
   out = zeros(Int,jsd,2)
   for i in 1:length(ns)
      out[ni[i,1]:ni[i,2],1] .= ns[i]
      out[ni[i,1]:ni[i,2],2] = collect(Int,-ns[i]:ns[i])
   end
   #[n k m]
   return out
end

struct psi
   F::Real #total angular momentum
   Ss::Array #spin angular momenta
   nfs::Array #symmetry fold of rotors
   ms::Array #mcalcs for respective rotors
   σs::Array #torsional symmetries for respective rotors
end
function qngen(ψ::psi)
   l = (2*ψ.f+1)*prod(2 .*ψ.Ss .+1)*prod(2 .*ψ.ms .+1)
   out = zeros(Int,l,length(ms)+length(Is)+2) #check this number of QNs
   out[:,1] = fill(Int(2*ψ.F),l)

end


function T10(sid::Int,p,ψb,ψk)
   #sid is the spin identity that for the T¹₀ operator to act on
   #p is the power that the operator will be raised to
   #ψb is the list of 
end


klist(N::Real) = collect(-N:N)
function klist(N::Array)
   out = zeros(Int(sum(2 .* N .+1)))
   si = 1
   fi = 0
   for n in N
      fi += Int(2*n+1)
      out[si:fi] = collect(-n:n)
      si += Int(2*n+1)
   end
   return out
end

function Δlist(J::Real,S::Real)
   return collect(abs(J-S):(J+S))
end
function Δlist(J::Array,S::Real)
   out = zeros(Int(length(J)*(2S+1)))
   sd = Int(2S+1)
   for i in eachindex(J)
      out[(i-1)*sd+1 : i*sd] = Δlist(J[i],S)
   end
   return out
end
function qngen(F,Ss)
   l = Int( (2*F+1)*prod(2 .*Ss .+1) )
   out = zeros(l, length(Ss) + 1)
   out[:,1] = fill(F,l)
   piece = F
   for i in 1:length(Ss)-1
      piece = Δlist(piece,Ss[i])
      out[:,i] = kron(piece,ones(Int16,Int(prod( 2 .*Ss[i+1:end] .+1)*(2F+1) )))
   end
   piece = Δlist(piece,Ss[end])
   out[:,end-1] = kron(piece, ones(Int(2F+1)))
   out[:,end] = klist(piece)
   return out
end

function Δlist16(dJ::Real,dS::Real)::Array{Int16} #takes doubled J and S!!!
   return collect(abs(dJ-dS):2:(dJ+dS))
end
function Δlist16_old(dJ::Array,dS::Real)::Array{Int16} #takes doubled J and S!!!
   sd = dS+1
   out = zeros(Int16,length(dJ)*sd)
   for i in eachindex(dJ)
      out[(i-1)*sd+1 : i*sd] = Δlist16(dJ[i],dS)
   end
   return out
end
function Δlist16(dJ::Array,dS::Real)::Array{Int16} #takes doubled J and S!!!
   out = reduce(vcat,Δlist16.(dJ,dS))
   return out
end

klist16(dN::Real)::Array{Int16} = collect(-dN:2:dN) #takes a doubled N!!!
function klist16(dN::Array)::Array{Int16} #takes a doubled N!!!
   out = zeros(sum(dN .+1))
   si = 1
   fi = 0
   for dn in dN
      fi += dn+1
      out[si:fi] = klist16(dn)
      si += dn+1
   end
   return out
end
d16(x::Real)::Int16 = Int16(2x)
degen(x::Real)::Int16 = Int16(2*x+1)
degen(x::Array)::Array{Int16} = Int16.(2 .* x .+1)
#This often works but needs a fix for F < sum(Ss)
function qngen16_prev(F::Real,Ss::Array)::Array{Int16,2}
   #if F ≤ sum(Ss)
   #   out = qngen16_safe(F,Ss)
   #else
      l = degen(F)*prod(degen(Ss))
      out = zeros(l, length(Ss) + 2)
      out[:,1] = fill(d16(F),l)
      piece = d16(F)
      for i in 1:length(Ss)-1
         piece = Δlist16(piece,d16(Ss[i]))
#         out[:,i+1] = kron(piece,ones( prod(degen(Ss[i+1:end]))*degen(F) ))
         @views out[:,i+1] = kron(piece,ones(Int16, Int16(l/length(piece)) ))
      end
      piece = Δlist16(piece,d16(Ss[end]))
      out[:,end-1] = kron(piece, ones(Int16, Int(l/length(piece))) )
      out[:,end] = klist16(piece)
      if 0 ∉ out[:,end]
         @warn "Non-integer N !!! Proceeding anyways but this is non-physical!"
      end
   #end
   return out
end
function qngen16(F::Real,Ss::Array)::Array{Int16,2}
   out = [d16(F)]
   for i in 1:length(Ss)
      piece = Δlist16(out[:,i],d16(Ss[i]))
      out = kron(out,ones(Int16, Int(length(piece)/length(out[:,i])) ))
      out = hcat(out, piece)
   end
   piece = klist16(out[:,end])
   out = kron(out, ones(Int16, Int(length(piece)/length(out[:,end])) ))
   out = hcat(out,piece)
   if 0 ∉ out[:,end]
      @warn "Non-integer N !!! Proceeding anyways but this is non-physical!"
   end
   return out
end

function qngen16_alt(F::Real,Ss::Array)::Array{Int16,2}
   out = Δlist16(d16(F),d16(Ss[1]))
   @show out
   for i in 2:length(Ss)
      piece = Δlist16(out[:,i-1],d16(Ss[i]))
      out = kron(out,ones(Int16, Int(length(piece)/length(out[:,i])) ))
      out = hcat(out, piece)
   end
   return out
end

