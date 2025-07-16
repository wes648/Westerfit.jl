
using LinearAlgebra, SparseArrays, BenchmarkTools
function ur(n::Int)::SparseMatrixCSC{Float64, Int}
   out = Diagonal(vcat(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(Diagonal(vcat(fill(√.5,n), 0.0, fill(√.5,n))))
   return sparse(out)
end
sgn(k::Real)::Int = iszero(k) ? 1 : sign(k)
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
δi(x::Real,y::Real)::Int = x==y

function pa_old(mc,a)
   ms = collect(-mc:mc)
   out = spdiagm(ms .^a)
   u = ur(mc)
   out = dropzeros!(u*out*u)
   return out
end
function pa_no(mc,a)
   ms = collect(-mc:mc)
   out = spdiagm(ms .^a)
   return out
end
function pa_new(mc,a)::SparseMatrixCSC{Float64, Int}
   out = map(x -> abs(x)^a, -mc:mc)
   if iseven(a)
      out = sparse(1:2mc+1, 1:2mc+1, out)
   else isodd(a)
      out = sparse(1:2mc+1, 2mc+1:-1:1, out)
   end
   return out#dropzeros!(out)
end

sand(a::AbstractArray,x::AbstractArray) = x' * a * x

function htorhc(mc,f,v)
   out = map(x -> f*x^2 + v, -mc:mc)
   out = sparse(reduce(vcat, [1:2mc+1; 1:2mc]),
                reduce(vcat, [1:2mc+1; 2:2mc+1]),
                vcat(out,fill(-0.5*v,2mc)))
   return out
end
function htorhc(nf,ms,f,v)
   out = map(x -> x^2 + v, -mc:mc)
   out = sparse(1:2mc+1, 1:2mc+1, out)
   out += coshc
   return out
end
function coshc(l,v,oddnf::Bool)
   if check
      out = spdiagm(1=>fill(0.5,l-1),-1=>fill(0.5,l-1))
   else
      out = spdiagm(2=>fill(0.5,l-2),-2=>fill(0.5,l-2))
   end
   u = ur(mc)
   out = sand(out,u)
   return out
end
function coshc(mc,v,a::Int)
   l = 2mc+1
   out = spdiagm(a=>fill(0.5,l-a),-a=>fill(0.5,l-a))
   u = ur(mc)
   out = dropzeros!(u*out*u)
   return out
end



function cos_old(mc,a)
   l = 2mc+1
   out = spdiagm(a=>fill(0.5,l-a),-a=>fill(0.5,l-a))
   u = ur(mc)
   out = dropzeros!(u*out*u)
   return out
end
function cos_dense(mc,a)
   l = 2mc+1
   out = diagm(a=>fill(0.5,l-a),-a=>fill(0.5,l-a))
   u = ud(mc)
   out = u*out*u
   return sparse!(out)
end


function cos_new(mc,a)
   l = 2mc+1
   c = spdiagm(a=>fill(0.5,l-a),-a=>fill(0.5,l-a))
   u = ur(mc)
   out = spzeros(l,l)
   out[1:mc,1:mc] = u[1:mc,:]*c*u[:,1:mc]
   out[mc+1:end,mc+1:end] = u[mc+1:end,:]*c*u[:,mc+1:end]
   return out
end
function cos_no(mc,a)
   l = 2mc+1
   out = spdiagm(a=>fill(0.5,l-a),-a=>fill(0.5,l-a))
   return out
end

function cos_b(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   out = spzeros(l,l)
   #term 1
   out .+= δi.((ms.+a)',ms)
   out .+= δi.((ms.-a)',ms)
   #term 2
   out .+= δi.((ms.+a)',-ms) .* sgn.(-mc:mc)
   out .+= δi.((ms.-a)',-ms) .* sgn.(-mc:mc)
   #term 3
   out .+= δi.((-ms.+a)',ms) .* sgn.(-mc:mc)'
   out .+= δi.((-ms.-a)',ms) .* sgn.(-mc:mc)'
   #term 4
   out .+= δi.((-ms.+a)',-ms) .* (sgn.(-mc:mc)' .*sgn.(-mc:mc))
   out .+= δi.((-ms.-a)',-ms) .* (sgn.(-mc:mc)' .*sgn.(-mc:mc))
   #clean up
   out .*= 0.25
   out[mc+1:end,mc+1] .*= √0.5
   out[mc+1,mc+1:end] .*= √0.5
   return dropzeros!(out)
end

function cos_b1(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   out = spzeros(l,l)
   #term 1
   out .+= δi.((ms.+a)',ms)
   out .+= δi.((ms.-a)',ms)
   #clean up
   out .*= 0.25
   out[mc+1:end,mc+1] .*= √0.5
   out[mc+1,mc+1:end] .*= √0.5
   return dropzeros!(out)
end
function cos_b2(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   out = spzeros(l,l)
   #term 2
   out .+= δi.((ms.+a)',-ms) .* sgn.(-mc:mc)
   out .+= δi.((ms.-a)',-ms) .* sgn.(-mc:mc)
   #clean up
   out .*= 0.25
   out[mc+1:end,mc+1] .*= √0.5
   out[mc+1,mc+1:end] .*= √0.5
   return dropzeros!(out)
end
function cos_b3(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   out = spzeros(l,l)
   #term 3
   out .+= δi.((-ms.+a)',ms) .* sgn.(-mc:mc)'
   out .+= δi.((-ms.-a)',ms) .* sgn.(-mc:mc)'
   #clean up
   out .*= 0.25
   out[mc+1:end,mc+1] .*= √0.5
   out[mc+1,mc+1:end] .*= √0.5
   return dropzeros!(out)
end
function cos_b4(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   out = spzeros(l,l)
   #term 4
   out .+= δi.((-ms.+a)',-ms) .* (sgn.(-mc:mc)' .*sgn.(-mc:mc))
   out .+= δi.((-ms.-a)',-ms) .* (sgn.(-mc:mc)' .*sgn.(-mc:mc))
   #clean up
   out .*= 0.25
   out[mc+1:end,mc+1] .*= √0.5
   out[mc+1,mc+1:end] .*= √0.5
   return dropzeros!(out)
end


function cos_w(mc,a)
   l = 2mc+1
   ms = abs.(-mc:mc)
   nm = collect(-mc:-1)
   pm = collect(0:mc)
   out = spzeros(l,l)
   #term 4
   out[1:mc,1:mc] .+= δi.((nm.+a)',nm) .* (sgn.(nm)' .*sgn.(nm))
   out[1:mc,1:mc] .+= δi.((nm.-a)',nm) .* (sgn.(nm)' .*sgn.(nm))

   out[mc+1:end,mc+1:end] .+= δi.((pm.+a)',pm) .* (sgn.(pm)' .*sgn.(pm))
   out[mc+1:end,mc+1:end] .+= δi.((pm.-a)',pm) .* (sgn.(pm)' .*sgn.(pm))

   #clean up
   out .*= 0.5
   out[mc+1+a,mc+1] .= √0.5
   out[mc+1,mc+1+a] .= √0.5
   return dropzeros!(out)
end
