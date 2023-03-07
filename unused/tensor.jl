
#=
include("./WIGXJPF.jl")
using LinearAlgebra, .WIGXJPF
include("./tensor.jl")
=#


function Δlist(J,S)
   max = Int(J+S)
   min = Int(abs(J-S))
   return collect(min:max)
end

function eh(x::Int)::Float64
   return √(x*(x+1))
end
#function fh(x::Int,y::Int)::Float64
#   out = □rt(x*(x+1) - y*(y+1))
#   return out
#end
function fh(x::Int,y::Int)::Float64
   out = □rt((x-y)*(x+y+1))
   return out
end
function gh(x,y)::Float64
   out = □rt((x-y)*(x-y-1.0))
   return out
end
function □rt(x)::Float64
   if x ≤ zero(x)
      return 0.0
   else
      return √x
   end
end
function θ(j,n,s)
   if s==0.5
      out = (n-j)/(j+0.5)
   else
      out = n*(n+1.0) + s*(s+1.0) - j*(j+1.0)
      out = out/(2.0*n*(n+1.0))
   end
   return out
end
function ϕ(j,n,s)
   if s==0.5
      out = 1.0
   else
      out = (n-j+s)*(n+j+s+1)*(s+j-n+1)*(n+j-s)
      out *= 1.0/((2.0*n-1.0)*(2.0*n+1.0))
      out = -√(out)/n
   end
   return out
end

δ(x::Int,y::Int)::Float64 = x==y
δ(x::Float64,y::Float64)::Float64 = x==y
δ(x::Float64,y::Int)::Float64 = x==convert(Float64,y)
δ(x::Int,y::Float64)::Float64 = convert(Float64,x)==y

T(l::Int,q::Int)::Int = l*(l+1) + q + 1

function kgen(n::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,1,2*n+1))
end
function kgen(n::Int,lb::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,1,lb))
end
function kgen(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = collect(-n:n)
      si = fi + 1
   end
   return kron(out,ones(Int,1,lb))
end
function ngen(n::Int)::Array{Int,2}
   return kron(fill(n,2*n+1),ones(Int,1,2*n+1))
end
function ngen(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = fill(n,2*n+1)
      si = fi + 1
   end
   return kron(out,ones(Int,1,lb))
end

nred(n::Int)::Float64 = √(n*(n+1)*(2*n+1))
nred(n::Float64)::Float64 = √(n*(n+1.0)*(2.0*n+1.0))
nnred(n::Int)::Float64 = n*(n+1)*(2*n+1)
nnred(n::Float64)::Float64 = n*(n+1.0)*(2.0*n+1.0)

function hrred(n::Int,l::Int)::Float64
   return wig6j(n,n,1,l,1,n)*nnred(n)*√(2*l+1)*(-1)^l
end
function hrelem(pr::Float64,l::Int,q::Int,nb,kb,nk,kk)::Array{Float64,2}
   @. return pr*wig3j(nb,l,nk,-kb,q,kk)*hrred(nk,l)*δ(nb,nk)*(-1)^(nk-kb)
end
function hrlpart(out,pr,l::Int,nb,kb,nk,kk)::Array{Float64,2}
   @simd for q in -l:l
      out .+= hrelem(pr[T(l,q)],l,q,nb,kb,nk,kk)
   end
   return out
end
function hrot(pr,n)
   nks = ngen(n)
   nbs = Matrix(transpose(nks))
   kks = kgen(n)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nbs,kbs,nks,kks)
   end
   return out
end
function hrot(pr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s,lb)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s,lb)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nbs,kbs,nks,kks)
   end
   return out
end

function nsred(l::Int,nb,nk)
   @. return 0.5*(nred(nk)*wig6j(1,1,l,nb,nk,nk) + nred(nb)*wig6j(1,1,l,nk,nb,nb))
end
function jsred(j,s,nb,nk)
   @. return wig6j(nk,s,j,s,nb,1)*jnred(nb,nk)*nnred(s)
end
function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)::Array{Float64,2}
   @. return pr*wig3j(nb,l,nk,-kb,q,kk)*√(2*l+1)*nsred(l,nb,nk)*jsred(j,s,nb,nk)*(-1)^(j+s-kb)
end
function srlpart(out,pr,l::Int,j,s,nb,kb,nk,kk)::Array{Float64,2}
   @simd for q in -l:l
      out .+= srelem(pr[T(l,q)],l,q,j,s,nb,kb,nk,kk)
   end
   return out
end
function hsr(pr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s,lb)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s,lb)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= srlpart(out,pr,l,j,s,nbs,kbs,nks,kks)
   end
   return out
end
function hrsr(pr,j,s)
   lb = convert(Int,(2.0*s+1.0)*(2.0*j+1.0))
   nks = ngen(j,s,lb)
   nbs = Matrix(transpose(nks))
   kks = kgen(j,s,lb)
   kbs = transpose(kks)
   out = zeros(size(kks))
   @simd for l in 0:2:2
      out .= hrlpart(out,pr,l,nbs,kbs,nks,kks)
      out .= srlpart(out,pr,l,j,s,nbs,kbs,nks,kks)
   end
   return out
end


function wigdiv(x::Array,s::Number)::Array
   if s==zeros(s)
      return 0.0 .* x
   else
      return x ./ wig3j(s,2,s,-s,0,1)
   end
end
function qured(j,s,nb,nk)
   @. return 0.25*jnred(nb,nk)*wig6j(j,s,nb,2,nk,s)
end
function quelem(pr,q,j,s,nb,kb,nk,kk)::Array{Float64,2}
   @. return pr*qured(j,s,nb,nk)*wig3j(nb,2,nk,-kb,q,kk)*(-1)^(nk+nb-kk+s+j+1)
end
function qulpart(out,pr,j,s,nb,kb,nk,kk)::Array{Float64,2}
   @simd for q in -l:l
      out .+= quelem(pr[T(l,q)],q,j,s,nb,kb,nk,kk)
   end
   out = wigdiv(out,s)
   return out
end

function jnred(j::Float64,n::Int)::Float64
   return √((2.0*j+1.0)*(2*n+1))
end
function jnred(j::Int,n::Int)::Float64
   return √((2*j+1)*(2*n+1))
end
function μred(s::Float64,jb::Float64,nb::Int,jk::Float64,nk::Int)::Float64
   return wig6j(nk,jk,s,jb,nb,1)*jnred(jb,nb)*jnred(jk,nk)*√(2*nb+1)
end
function μelem(pr::Float64,q,s::Float64,jb::Float64,nb,kb,jk::Float64,nk,kk)::Array{Float64,2}
   @. return pr*wig3j(nb,1,nk,-kb,q,kk)*μred(s,jb,nb,jk,nk)
end
function intmat(μs,s,jb,jk)
   lb = Int((2.0*s+1.0)*(2.0*jb+1.0))
   lk = Int((2.0*s+1.0)*(2.0*jk+1.0))
   nks = transpose(ngen(jk,s,lb))
   kks = transpose(kgen(jk,s,lb))
   nbs = ngen(jb,s,lk)
   kbs = kgen(jb,s,lk)
   out = zeros(lb,lk)
   @simd for q in -1:1
      out .+= μelem(μs[T(1,q)],q,s,jb,nbs,kbs,jk,nks,kks)
   end
   return out
end


A = 3.0
B = 1.5
C = 1.0
Dab = 0.1

prm = [-(A+B+C)/√3; 0; 0; 0; (B-C)/2; Dab; (2*A-B-C)/√6; -Dab; (B-C)/2]



function cart2sphr(inp::Array{Float64,2})::Array{Float64,1}
   out = zeros(9)
   out[1] = -sum(diag(inp))/√3
   out[2] = 0.5*(inp[1,2] - inp[2,1])
   out[3] = 0.0
   out[4] = 0.5*(inp[1,2] - inp[2,1])
   out[5] = 0.5*(inp[2,2] - inp[3,3])
   out[6] = 0.5*(inp[1,2] + inp[2,1])
   out[7] = (3.0*inp[1,1] - sum(diag(inp)))/√6
   out[8] = -0.5*(inp[1,2] + inp[2,1])
   out[9] = 0.5*(inp[2,2] - inp[3,3])
   return out
end
function cart2sphr(inp::Array{Float64,1})::Array{Float64,1}
   out = zeros(3)
   out[1] = inp[2]/√2
   out[2] = inp[1]
   out[3] = -inp[2]/√2
   return out
end


function nzop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. δ(nb,nk)*δ(kb,kk)*kk^p
end
function ntop(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. δ(nb,nk)*δ(kb,kk)*eh(nk)^p
end
function npel(nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. δ(nb,nk)*δ(kb-1,kk)*fh(nk,kk)
end
function nmel(nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. δ(nb,nk)*δ(kb+1,kk)*fh(nk,kk-1)
end
function npmp(p::Int,nb::Array{Int,2},kb::Array{Int,2},
              nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return npel(nb,kb,nk,kk)^p .+ nmel(nb,kb,nk,kk)^p
end
function rotop(pr::Float64,t::Int,z::Int,x::Int,
               nb::Array{Int,2},kb::Array{Int,2},
               nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   tz = ntop(t,nb,kb,nk,kk).*nzop(z,nb,kb,nk,kk) #./ 2.0
   pm = npmp(x,nb,kb,nk,kk)
   return pr .* (pm*tz .+ tz*pm) ./ 4.0
end
function nsop(p::Int,j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return @. 0.5*δ(nb,nk)*δ(kb,kk)*(eh(j) - eh(nk) - eh(s))^p
end
function szen(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb,nk)*δ(kb,kk)*k*θ(j,nk,s)
end
function szem(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb-1,nk)*δ(kb,kk)*√(nk^2 - kk^2)*ϕ(j,nk,s)
end
function szep(j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})
   return @. δ(nb+1,nk)*δ(kb,kk)*√(nb^2 - kk^2)*ϕ(j,nb,s)
end
function szop(p::Int,j::Float64,s::Float64,nb::Array{Int,2},
              kb::Array{Int,2},nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   return (szen(j,s,nb,kb,nk,kk) .+ szem(j,s,nb,kb,nk,kk) .+ szep(j,s,nb,kb,nk,kk))^p
end
function rsrop(pr::Float64,a::Int,b::Int,c::Int,d::Int,e::Int,
               j,s,nb::Array{Int,2},kb::Array{Int,2},
               nk::Array{Int,2},kk::Array{Int,2})::Array{Float64,2}
   tz = ntop(a,nb,kb,nk,kk).*nzop(b,nb,kb,nk,kk).*nsop(d,j,s,nb,kb,nk,kk) #these all commute!
   pm = npmp(c,nb,kb,nk,kk)
   sz = szop(e,j,s,nb,kb,nk,kk)
   return pr .* (tz*sz*pm + pm*sz*tz) ./ 4.0
end


nk = ngen(10.,0.,21)
kk = kgen(10.,0.,21)
nb = Matrix(transpose(nk))
kb = Matrix(transpose(kk))
@time rotop(1.0,2,2,2,nb,kb,nk,kk)