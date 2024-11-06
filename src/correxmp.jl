

using Random, Distributions
using Plots, LaTeXStrings

Plots.default(fontfamily = "Computer Modern")

include("./main.jl")

function hr(n)
   out = ones(2*n+1,2*n+1) - diagm(1=>fill(1,2*n),-1=>fill(1,2*n))
   out .*= hrtest(n)
   return Symmetric(Matrix(out))
end

function hr(n,κ)
   a = 3e3
   b = (κ+2)*1e3
   c = 1e3   
   if κ < -0.0001 #prolate Ir
      bk = a - 0.5*(b+c)
      bn = 0.5*(b+c)
      bp = 0.25*(b-c)
   elseif κ > 0.0001 #oblate IIIr
      bk = c - 0.5*(a+b)
      bn = 0.5*(a+b)
      bp = 0.25*(a-b)
   else # accidental symmetric IIr
      bk = b - 0.5*(c+a)
      bn = 0.5*(c+a)
      bp = 0.25*(c-a)
   end
   out = bk * nz2(n)
   out += bn * n2(n)
   out += bp * npm2(n)
   return Symmetric(Matrix(out))
end
rms(x) = BLAS.nrm2(x)/√length(x)


function quickeigs(nmax,κ)
   vals = zeros(sum(2 .*collect(0:nmax) .+ 1))
   vecs = zeros(2*nmax+1,length(vals))
   si = 2
   fi = 4
   vals[1] = eigvals(hr(0,κ))[1]
   vecs[1,1] = eigvecs(hr(0,κ))[1]
   for n in 1:nmax
      vals[si:fi],vecs[1:(2*n+1),si:fi] = eigen(hr(n,κ))
      si += 2*n+1
      fi += 2*(n+1)+1
   end
   return vals, vecs
end

function n2(n)
   diagm(fill(n*(n+1),2*n+1 ))
end
function nz2(n)
   diagm(collect(-n:n).^2)
end
function npm2(n)
   nk = ngen(n)
   nb = permutedims(nk)
   kk = kgen(n)
   kb = permutedims(kk)
   return npmp(2,nb,kb,nk,kk)
end

function t400(n)
   diagm(fill((n*(n+1))^2, 2*n+1))
end
function t040(n)
   diagm(collect(-n:n).^4)
end
function t220(n)
   diagm(fill(n*(n+1), 2*n+1))*diagm(collect(-n:n).^2)
end
function t202(n)
   nk = ngen(n)
   nb = permutedims(nk)
   kk = kgen(n)
   kb = permutedims(kk)
   p1 = npmp(2,nb,kb,nk,kk)
   p2 = n2(n)
   return p1*p2 + p2*p1
end
function t022(n)
   nk = ngen(n)
   nb = permutedims(nk)
   kk = kgen(n)
   kb = permutedims(kk)
   p1 = npmp(2,nb,kb,nk,kk)
   p2 = nz2(n)
   return p1*p2 + p2*p1
end
function t004(n)
   nk = ngen(n)
   nb = permutedims(nk)
   kk = kgen(n)
   kb = permutedims(kk)
   return npmp(4,nb,kb,nk,kk)
end

function bkder(nmax,vecs)
   der = zeros(sum(2 .*collect(0:nmax) .+ 1))
   der[1] = 0.0
   si = 2
   fi = 4
   for n in 1:nmax
      v = vecs[1:(2*n+1),si:fi]
      der[si:fi] = diag(v' * nz2(n) * v)
      si += 2*n+1
      fi += 2*(n+1)+1
   end
   return der
end
function bnder(nmax,vecs)
   der = zeros(sum(2 .*collect(0:nmax) .+ 1))
   der[1] = 0.0
   si = 2
   fi = 4
   for n in 1:nmax
      der[si:fi] = diag(vecs[1:(2*n+1),si:fi]' * n2(n) * vecs[1:(2*n+1),si:fi])
      si += 2*n+1
      fi += 2*(n+1)+1
   end
   return der
end
function bpder(nmax,vecs)
   der = zeros(sum(2 .*collect(0:nmax) .+ 1))
   der[1] = 0.0
   si = 2
   fi = 4
   for n in 1:nmax
      der[si:fi] = diag(vecs[1:(2*n+1),si:fi]' * npm2(n) * vecs[1:(2*n+1),si:fi])
      si += 2*n+1
      fi += 2*(n+1)+1
   end
   return der
end

function μgen(n)
   nk = ngen(n)
   nb = permutedims(nk)
   kk = kgen(n)
   kb = permutedims(kk)
   out = zeros(2*n+1,2*n+1)
   @. out = wig3j(nb,1,nk,-kb,-1,kk) + wig3j(nb,1,nk,-kb,0,kk) + wig3j(nb,1,nk,-kb,1,kk)
   return sparse!(out)
end
function μgen(k,b)
   nk = ngeni(k,2*b+1)
   nb = permutedims(ngen(b,2*k+1))
   kk = kgen(k,2*b+1)
   kb = permutedims(kgen(b,2*k+1))
   out = zeros(2*k+1,2*b+1)
   @. out = wig3j(nb,1,nk,-kb,-1,kk) + wig3j(nb,1,nk,-kb,0,kk) + wig3j(nb,1,nk,-kb,1,kk)
   return sparse!(out)
end
function allowints(nmax,vecs)
   mat = spzeros(size(vecs,2),size(vecs,2))
   for nk in 0:nmax
   sk,fk = jinds(nk,0)
   vk = vecs[1:(2*nk+1),sk:fk]
   for nb in 0:nmax
      sb,fb = jinds(nb,0)
      vb = vecs[1:(2*nb+1),sb:fb]
      part = vb' * μmat(ones(3),0.0,Float64(nb),Float64(nk))
      part *= vk
      mat[sb:fb,sk:fk] .= part
   end
   end
   return droptol!(mat,1e-9)
end
function allowtras(nmax::Real,vals::Array{Float64},vecs::Array{Float64,2})
   out = spzeros(length(vals),length(vals))
   ro,co,v = findnz(allowints(nmax,vecs))
   for i in 1:length(ro)
      r = ro[i]
      c = co[i]
      out[r,c] = vals[r] - vals[c]
   end
   out .*= (out .> 0.0)
   return droptol!(sparse!(out), 1e-9)
end
function allowtras(nmax::Real,κ::Real)
   vals,vecs = quickeigs(nmax,κ)
   return allowtras(nmax,vals,vecs)
end

function gradf(f::Function,nmax,vecs)
   out = zeros(size(vecs,2))
   for n in 1:nmax
      si,fi = jinds(n,0)
      vs = vecs[1:(2*n+1), si:fi]
      out[si:fi] .= diag(vs' * f(n) *vs) 
   end
   return out
end
function allgrads(fs::Array{Function},nmax,vecs)
   out = zeros(size(vecs,2),length(fs))
   for i in eachindex(fs)
      out[:,i] = gradf(fs[i],nmax,vecs)
   end
   return out
end
function jacobians(fs::Array{Function},nmax,rows,cols,grads)
   out = zeros(length(rows),length(fs))
   for f in eachindex(fs)
   for i in eachindex(rows)
      r = rows[i]
      c = cols[i]
      out[i,f] = grads[r,f] - grads[c,f]
   end
   end
   return out
end

J(nmax,vs) = [bkder(nmax,vs) bnder(nmax,vs) bpder(nmax,vs)]
H(nmax,vs) = J(nmax,vs)' * J(nmax,vs)
H(j) = j' * j

function hess2corr(h)
   corr = zero(h)
   for i in 1:size(h,1)
   for j in 1:size(h,2)
      corr[i,j] = h[i,j] / √(h[i,i]*h[j,j])
   end
   end
   return corr
end

function corrgen(flst,nmax,κ)
   val,vec = quickeigs(nmax,κ)
   row,col,v = findnz(allowtras(nmax,val,vec))
   grad = allgrads(flst,nmax,vec)
   jac = jacobians(flst, nmax,row,col,grad)
   return hess2corr(H(jac))
end

function bigtest(nmax,s)
   fs = [nz2;n2;npm2;t400;t220;t040;t202;t022;t004]
   κs = rand(Uniform(-1,1),s)
   out = zeros(9,9,s)
   #out[:,1] = κs
   @threads for i in eachindex(κs)
      out[:,:,i] = corrgen(fs,nmax,κs[i])
   end
   return out, κs
end

function bkplot(nmax,s)
   res, κ = bigtest(nmax,s)
   p1 = histogram(κ,bins=15)
   p2 = scatter(κ,res[2,1,:],title=L"B_{N}")
   p3 = scatter(κ,res[3,1,:],title=L"B_{\pm}")
   p4 = scatter(κ,res[4,1,:],title=L"T_{400}")
   p5 = scatter(κ,res[5,1,:],title=L"T_{220}")
   p6 = scatter(κ,res[6,1,:],title=L"T_{040}")
   p7 = scatter(κ,res[7,1,:],title=L"T_{202}")
   p8 = scatter(κ,res[8,1,:],title=L"T_{022}")
   p9 = scatter(κ,res[9,1,:],title=L"T_{004}")
   plot(p1,p2,p3,p4,p5,p6,p7,p8,p9, layout=(3,3), legend=false, dpi=400,xlab=L"\kappa",
      ylab="Corr "*L"B_{K}")
end   
function bnplot(nmax,s)
   res, κ = bigtest(nmax,s)
   p1 = histogram(κ,bins=15)
   p2 = scatter(κ,res[1,2,:],title=L"B_{K}")
   p3 = scatter(κ,res[3,2,:],title=L"B_{\pm}")
   p4 = scatter(κ,res[4,2,:],title=L"T_{400}")
   p5 = scatter(κ,res[5,2,:],title=L"T_{220}")
   p6 = scatter(κ,res[6,2,:],title=L"T_{040}")
   p7 = scatter(κ,res[7,2,:],title=L"T_{202}")
   p8 = scatter(κ,res[8,2,:],title=L"T_{022}")
   p9 = scatter(κ,res[9,2,:],title=L"T_{004}")
   plot(p1,p2,p3,p4,p5,p6,p7,p8,p9, layout=(3,3), legend=false, dpi=400,xlab=L"\kappa",
      ylab="Corr "*L"B_{N}")
end  
function bpplot(nmax,s)
   res, κ = bigtest(nmax,s)
   p1 = histogram(κ,bins=15)
   p2 = scatter(κ,res[1,3,:],title=L"B_{K}")
   p3 = scatter(κ,res[2,3,:],title=L"B_{N}")
   p4 = scatter(κ,res[4,3,:],title=L"T_{400}")
   p5 = scatter(κ,res[5,3,:],title=L"T_{220}")
   p6 = scatter(κ,res[6,3,:],title=L"T_{040}")
   p7 = scatter(κ,res[7,3,:],title=L"T_{202}")
   p8 = scatter(κ,res[8,3,:],title=L"T_{022}")
   p9 = scatter(κ,res[9,3,:],title=L"T_{004}")
   plot(p1,p2,p3,p4,p5,p6,p7,p8,p9, layout=(3,3), legend=false, dpi=400,xlab=L"\kappa",
      ylab="Corr "*L"B_{\pm}")
end  