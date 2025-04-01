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
   out = Diagonal(ms .^a)
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

function pafull(mk,a)
   mb = mk'
   out = spzeros(length(mk),length(mk))
   out += @. mk^a *δi(mb,mk)
   out -= @. sgn(mk)*(-mk)^a *δi(mb,-mk)
   out -= @. sgn(mb)*mk^a *δi(-mb,mk)
   out += @. sgn(mk)*sgn(mb)*(-mk)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end

function pa_test(mk,a)
   mb = mk'
   ms = mk#abs.(mk)
   out = spzeros(length(mk),length(mk))
   out += @. ms^a *δi(mb,mk)
   out += @. sgn(mk)* (-ms)^a *δi(-mb,mk)
   out += @. sgn(mb)* ( ms)^a *δi(mb,-mk)
   out += @. sgn(mk)*sgn(mb)*(-ms)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end
function pa_test14(mk,a)
   mb = mk'
   ms = abs.(mk)
   out = spzeros(length(mk),length(mk))
   out += @. ms^a *δi(mb,mk)
   #out += @. sgn(mk)* (-ms)^a *δi(mb,-mk)
   #out += @. sgn(mb)* ms^a *δi(-mb,mk)
   out += @. sgn(mk)*sgn(mb)*(-ms)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end
function pa_test23(mk,a)
   mb = mk'
   ms = abs.(mk)
   out = spzeros(length(mk),length(mk))
   #out += @. ms^a *δi(mb,mk)
   out += @. sgn(mk)* (-ms)^a *δi(mb,-mk)
   out += @. sgn(mb)* ms^a *δi(-mb,mk)
   #out += @. sgn(mk)*sgn(mb)*(-ms)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end
function pa_test2(mk,a)
   mb = mk'
   ms = abs.(mk)
   out = spzeros(length(mk),length(mk))
   #out += @. ms^a *δi(mb,mk)
   out += @. sgn(mk)* (-ms)^a *δi(mb,-mk)
   #out += @. sgn(mb)* ms^a *δi(-mb,mk)
   #out += @. sgn(mk)*sgn(mb)*(-ms)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end
function pa_test3(mk,a)
   mb = mk'
   ms = abs.(mk)
   out = spzeros(length(mk),length(mk))
   #out += @. ms^a *δi(mb,mk)
   #out += @. sgn(mk)* (-ms)^a *δi(mb,-mk)
   out += @. sgn(mb)* ms^a *δi(-mb,mk)
   #out += @. sgn(mk)*sgn(mb)*(-ms)^a *δi(-mb,-mk)
   out = 0.5 .* sparse(out)
   return out
end


function pa_old(mc,a)
   ms = collect(-mc:mc)
   out = Diagonal(ms .^a)
   u = ur(mc)
   out = dropzeros!(u*out*u)
   return out
end

function pa_new(mc,a)::SparseMatrixCSC{Float64, Int}
   ms = collect(-mc:mc)
   if iseven(a)
      out = spdiagm(ms .^ a)
   else
      out = rotl90(spdiagm(abs.(ms).^a))
   end
   return dropzeros!(out)
end