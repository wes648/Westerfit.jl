
function θ(j,n,s)
   if s==0.5
      out = (n-j)/(j+0.5)
   #elseif n==0.0
   #   out = 0.0
   else
      out = n*(n+1.0) + s*(s+1.0) - j*(j+1.0)
      out = out/(2.0*n*(n+1.0))
   end
   return out
end
function ϕ(j,n,s)
   if s==0.5
      out = -1.0/(j+0.5)
   else
      out = (n-j+s)*(n+j+s+1)*(s+j-n+1)*(n+j-s)
      out *= 1.0/((2.0*n-1.0)*(2.0*n+1.0))
      out = -sqrt(out)/n
   end
   return out
end

function fh(x,y)::Float64
   out = sqrt(x*(x+1.0) - y*(y+1.0))
   return out
end
function gh(x,y)::Float64
   out = sqrt((x-y)*(x-y-1.0))
   return out
end

function Ediag(pr,J,S,N,K,m,σ)
   out = Htor0m(pr,N,K,m,σ)
   out += Hr0K(pr,N,K)
   out += Hs0K0N(pr,J,N,S,θ(J,N,S),K)
   return out
end

function Htor0m(pr,N,K,m,σ)
   #<m K|Htor|m K> =
   out = @. pr[5]*m^2 + 2*pr[6]*m*K + pr[7]*0.5
end
function Htor1m(pr,N,K,m,σ)
   #<m+1 K|Htor|m K> = -V3/4
   out = @. -pr[7]*0.25 + 0.0*m*K
end
function Hr0K(pr,N,K)
   #<N K |Hrot| N K>
   #   out = @. (0.5*(pr[2] + pr[3])+ pr[13]*N*(N+1.0))*N*(N+1.0) + (pr[1] - 0.5*(pr[2] + pr[3]))*K^2
   out = @. (pr[2]+ pr[13]*N*(N+1.0))*N*(N+1.0) + pr[1]*K^2
end
function Hr1K(pr,N,K)
   #<N K+1 |Hrot|N K>
   out = @. pr[4]*sqrt(N*(N+1.0)-K*(K-1.0))*(K-0.5)
end
function Hr2K(pr,N,K)
   #<N K+2 |Hrot| N K>
   out = @. pr[3]*sqrt((N*(N+1.0) - K*(K - 1.0))*(N*(N+1.0) - (K - 1.0)*(K - 2.0)))
end
function Hrot(pr,N)
   if N == zero(N)
      return spzeros(Float64,1,1)
   else
      karray = collect(Float64,-N:N)
      nd = length(karray)
      ondiags = Hr0K(pr,N,karray)
      of1diag = Hr1K(pr,N,karray[2:end])
      of2diag = Hr2K(pr,N,karray[3:end])
      out = spdiagm(nd,nd,0=>ondiags,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
      return out
   end
end

function Htor(pr,mcalc,N,σ)
   tnp = convert(Int,2*N+1)
   if mcalc == zero(mcalc)
      return spzeros(tnp,tnp)
   else
   tmp = convert(Int,2*mcalc+1)
   ks = kron(ones(Float64,tmp),collect(Float64,-N:N))
   ms = kron(NFOLD .* collect(Float64,-mcalc:mcalc) .+ σ, ones(Float64,tnp))
   ondiags = Htor0m(pr,N,ks,ms,σ)
   of1diags = Htor1m(pr,N,ks[tnp+1:end],ms[tnp+1:end],σ)
   out = spdiagm(tnp*tmp,tnp*tmp,0=>ondiags,tnp=>of1diags,-tnp=>of1diags)
   return out
   end
end
function Htor(pr,mcalc,j,s,σ)
   tnp = convert(Int,(2*j+1)*(2*s+1))
   if mcalc == zero(mcalc)
      return spzeros(tnp,tnp)
   else
   tmp = convert(Int,2*mcalc+1)
   ks = zeros(0)
   for n in Δlist(j,s)
      ks = vcat(ks,collect(Float64,-n:n))
   end
   ks = kron(ones(Float64,tmp),ks)
   ms = kron(NFOLD .* collect(Float64,-mcalc:mcalc) .+ σ, ones(Float64,tnp))
   ondiags = Htor0m(pr,j,ks,ms,σ)
   of1diags = Htor1m(pr,j,ks[tnp+1:end],ms[tnp+1:end],σ)
   out = spdiagm(tnp*tmp,tnp*tmp,0=>ondiags,tnp=>of1diags,-tnp=>of1diags)
   return out
   end
end

function Htr(pr,N,m,σ)
   nd = 2*N+1
   md = 2*m+1
   out = spzeros(Float64,md*nd,md*nd)
   out .= kron(eye(md),Hrot(pr,N))
   out += Htor(pr,m,N,σ)
   return out
end
function Hs0K0N(pr,J,N,S,thet,K)
#   <J S N K| Hsr |J S N K>
   out  = @. -0.5*pr[8]*(J*(J+1.0)-N*(N+1)-S*(S+1.0)) + thet*(pr[9]*(3.0*K^2 - N*(N+1)))
end
function Hs1K0N(pr,J,N,thet,K)
#   <J S N K-1| Hsr |J S N K>
   out = @. pr[11]*(K-0.5)*fh(N,K-1.0)*thet
end
function Hs2K0N(pr,J,N,thet,K)
   #   <J S N K-2| Hsr |J S N K>
   out = @. thet*0.5*pr[10]*sqrt((N*(N+1)-K*(K-1.0))*(N*(N+1)-(K-1.0)*(K-2.0)))
end
function Hsm2K1N(pr,N::Float64,K::Array{Float64,1},ϕ)#::Array{Float64,1}
#   <J S N+1 K-2| Hsr |J S N K>
   #out = @. -0.25*(srp[3]+srp[12]*(K*(N-K)+(K-2)*(N-K+2)))*f(N,K-1)*g(N,-K+1.0)
   out = @. 0.25*pr[10]*fh(N+1.0,K-2.0)*gh(N+1.0,K-1.0).*ϕ
end
function Hsm1K1N(pr,N::Float64,K::Array{Float64,1},ϕ)#::Array{Float64,1}
#   <J S N+1 K-1| Hsr |J S N K>
   out = @. 0.25*pr[11]*(N+2.0*K)*gh(N+1.0,K-1.0).*ϕ
end
function Hs0K1N(pr,N::Float64,K::Array{Float64,1},ϕ)#::Array{Float64,1}
#   <J S N+1 K| Hsr |J S N K>
   out = @. 1.5*pr[9]*K*sqrt((N+1.0)^2 - K^2).*ϕ
end
function Hs1K1N(pr,N::Float64,K::Array{Float64,1},ϕ)#::Array{Float64,1}
#   <J S N+1 K+1| Hsr |J S N K>
   out = @. 0.25*pr[11]*(N-2.0*K)*gh(N+1.0,-K-1.0).*ϕ
end
function Hs2K1N(pr,N::Float64,K::Array{Float64,1},ϕ)#::Array{Float64,1}
#   <J S N+1 K+2| Hsr |J S N K>
   out = @. -0.25*pr[10]*fh(N+1.0,K+1.0)*gh(N+1.0,-K-1.0).*ϕ
end
#quick hyperfine implementation so I can stop thinking about it
function hqqelem(χ,f,i,jb,j,k,q)
   kb = q+k
   out = wig3j(jb,2,jk,-kb,q,k)*Tχ(q)
   if out == zero(out)
      return out
   else
      out *= 0.25*χ[abs(q)]*√((2.0*jb+1.0)*(2.0*j+1.0))*(-1)^(jb+jk-kb+i+f+1)
   end
end
function hqqmat(pr,f,i,jb,j)
   χ = pr[13:15]
   out = spzeros(float64,Int(2*jb)+1,Int(2*j)+1)
   fac = wig6j(f,i,j,2,jb,i)
   if (fac == zero(fac))||(i==zero(i))
      return out
   else
      fac /= wig3j(i,2,i,-i,0,i)
      Δ = Int(jb-j)
      karray = collect(Float64,-j:j)
      pm2 = hqqelem(χ,f,i,jb,j,karray[2:end],-2)
      pm1 = hqqelem(χ,f,i,jb,j,karray,-1)
      p0  = hqqelem(χ,f,i,jb,j,karray,0)
      p1  = hqqelem(χ,f,i,jb,j,karray,1)
      p2  = hqqelem(χ,f,i,jb,j,karray[1:end-1],2)
      out .= spdiagm(-2+Δ=>pm2,-1+Δ=>pm1,0+Δ=>p1,0+Δ=>p1,2+Δ=>p2,)
      return out
   end
end

#Matrix Builder
function Hspi0N(pr,J,S,N)
   if N == zero(N)
      #ondiags = Hs0K0N(J,N,S,Nt2,thet,[0.0],srprms)
      out = spzeros(Float64,1,1)#diagm(0=>ondiags)
   else
      thet = θ(J,N,S)
      Nt2 = N*(N+1.0)
      karray = collect(Float64,-N:N)
      ondiags = Hs0K0N(pr,J,N,S,thet,karray)
      of1diag = Hs1K0N(pr,J,N,thet,karray[2:end])
      of2diag = Hs2K0N(pr,J,N,thet,karray[3:end])
      out = spdiagm(length(karray),length(karray),0=>ondiags,1=>of1diag,2=>of2diag,-1=>of1diag,-2=>of2diag)
   end
   return out
end
function Hspi1N(pr,J::Float64,S::Float64,Nl)
   jϕ = ϕ(J,Nl+1.0,S)
   karray = collect(Float64,-Nl:Nl)
   pm2 = Hsm2K1N(pr,Nl,karray[2:end],jϕ)
   pm1 = Hsm1K1N(pr,Nl,karray,jϕ)
   p0 = Hs0K1N(pr,Nl,karray,jϕ)
   p1 = Hs1K1N(pr,Nl,karray,jϕ)
   p2 = Hs2K1N(pr,Nl,karray[1:end-1],jϕ)
   Spimat = spdiagm(2*Int(Nl)+1,2*Int(Nl)+3, -1=>pm2, 0=>pm1, 1=>p0, 2=>p1, 3=>p2)
   return Spimat
end

function Htsr0N(pr,j,s,n,mcalc,σ)
   marray = NFOLD.* collect(Float64,-mcalc:mcalc) .+ σ
   karray = collect(Float64,-n:n)
   if n == zero(n)
      mat = spzeros(Float64,length(marray),length(marray))#diagm(0=>ondiags)
   else
      mat = pr[12]*θ(j,n,s)*0.0
      mat = mat .* karray
      mat = kron(marray,mat)
      mat = spdiagm(0=>mat)
   end
   return mat
end
function Htsr1N(pr,j::Float64,s::Float64,nl,mcalc,σ)
   karray = collect(Float64,-nl:nl)
   marray = NFOLD .* collect(Float64,-mcalc:mcalc) .+ σ
   mat = spzeros(Float64,0,0)
   for i in 1:length(marray)
      p1 = ϕ(j,nl+1.0,s)*pr[12]*marray[i]*0.0
      p1 = p1 .* sqrt.( (nl+1.0)^2 .+ karray .^2)
      part = spdiagm((2*Int(nl)+1),(2*Int(nl)+3), 1=>p1)
      mat = cat(mat,part,dims=(1,2))
   end
   return mat
end

function Hsr(pr,J,S)
   ns, nd, ni, jd = srprep(J,S)
   out = spzeros(Float64,jd,jd)
   out[1:nd[1],1:nd[1]] = Hrot(pr,ns[1]) + Hspi0N(pr,J,S,ns[1])
   for i in 2:length(ns)
      n = ns[i]
      n1part = Hspi1N(pr,J,S,n-1.0)
      out[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]] = n1part
      out[   ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]] = Hrot(pr,n) + Hspi0N(pr,J,S,n)
      out[   ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(n1part)
   end
   return out
end
function Htsr(pr,J,S,mcalc,σ)#OLD
   md = 2*mcalc + 1
   ns, nd, ni, jd = srprep(J,S,md)
#   ni = ni
   ni[1,1] = 1 #investigate why this is needed should be definitional to srprep but isn't???
   out = spzeros(Float64,md*jd,md*jd)
   out[1:ni[1,2],1:ni[1,2]] = kron(eye(md),Hrot(pr,ns[1]) + Hspi0N(pr,J,S,ns[1]))
   out[1:ni[1,2],1:ni[1,2]] += Htor(pr,mcalc,ns[1],σ) + Htsr0N(pr,J,S,ns[1],mcalc,σ)
   for i in 2:length(ns)
   n = ns[i]
   n1part = kron(eye(md),Hspi1N(pr,J,S,n-1.0)) + Htsr1N(pr,J,S,n-1.0,mcalc,σ)
   @inbounds out[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]] = n1part
   @inbounds out[   ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]] = kron(eye(md),
               Hrot(pr,n) + Hspi0N(pr,J,S,n)) + Htor(pr,mcalc,n,σ) + Htsr0N(pr,J,S,n,mcalc,σ)
   @inbounds out[   ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(n1part)
   end
   return out
end

function Htsr0Nv(pr,j,s,n)
   karray = collect(Float64,-n:n)
   if n == zero(n)
      mat = spzeros(Float64,1,1)
   else
      mat = pr[12]*θ(j,n,s)
      mat = mat .* karray
      mat = spdiagm(0=>mat)
   end
   return mat
end
function Htsr1Nv(pr,j,s,nl)
   karray = collect(Float64,-nl:nl)
   p1 = ϕ(j,nl+1.0,s)*pr[12]
   p1 = p1 .* sqrt.( (nl+1.0)^2 .+ karray .^2)
   mat = spdiagm((2*Int(nl)+1),(2*Int(nl)+3), 1=>p1)
   return mat
end

function Htsrmat2(pr,j,s,mcalc,σ)
   ns, nd, ni, jd = srprep(j,s)
   srpart = spzeros(Float64,jd,jd)
   tspart = spzeros(Float64,jd,jd)
   srpart[1:nd[1],1:nd[1]] = Hrot(pr,ns[1]) + Hspi0N(pr,j,s,ns[1])
   tspart[1:nd[1],1:nd[1]] = Htsr0Nv(pr,j,s,ns[1])
   trpart = Htor(pr,mcalc,j,s,σ)
   for i in 2:length(ns)
      n = ns[i]
      srn1part = Hspi1N(pr,j,s,n-1.0)
   @inbounds srpart[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]] = srn1part
   @inbounds srpart[   ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]] = Hrot(pr,n)+ Hspi0N(pr,j,s,n)
   @inbounds srpart[   ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(srn1part)
      n1part = Htsr1Nv(pr,j,s,n-1.0)
   @inbounds tspart[ni[i-1,1]:ni[i-1,2],   ni[i,1]:ni[i,2]] = n1part
   @inbounds tspart[   ni[i,1]:ni[i,2],   ni[i,1]:ni[i,2]] = Htsr0Nv(pr,j,s,n)
   @inbounds tspart[   ni[i,1]:ni[i,2],ni[i-1,1]:ni[i-1,2]] = transpose(n1part)
   end
   marray = NFOLD.* collect(Float64,-mcalc:mcalc) .+ σ
   tspart = kron(diagm(0=>marray),tspart)
   out = kron(eye(2*mcalc+1),srpart) + tspart + trpart
   return out
end
