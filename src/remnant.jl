"""
This is the remnant code from userdefop. It needs to either be moved elsewhere
or fully replaced
"""
function indpuller(vtm,mc,σt,jsd)
   a = vt2m(vtm,σt)
   b = vt2m(vtm-1,σt)
   ref = sort([a; b])
   if σt != 2
      c = jsd*(mc + ref[1]) + 1
      d = jsd*(mc + ref[2] +1)
   else
      c = jsd*(mc + ref[1] - (vtm>zero(vtm))) + 1
      d = jsd*(mc + ref[2] + 1 - (vtm>zero(vtm)))
   end
   return collect(c:d)
end
function jvdest(j,s,vtm)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd,fnd
end
function jvlinds(j,s,vtm)
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return collect(snd:fnd)
end
function jinds(j,s,m,σt)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*
                       sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*
                       sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd,fnd
end
function jlinds(j,s,m,σt)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*
                       sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (2*m+1+δ(σt,2))*(2*s+1)*
                       sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return collect(snd:fnd)
end

function qnlab(n,nf,m,σ)
   nd = Int(2*n+1)
   md = Int(2*m+(σtype(nf,σ)==2)+1)
   narray = fill(n,nd*md)
   karray = kron(ones(Int,md),collect(Int,-n:n))
   kcarray = k2kc.(narray,karray)
   marray = kron(msbuilder(nf,m,σ),ones(Int,nd))
   σarray = fill(σ,nd*md)
   out = hcat(narray,abs.(karray),kcarray,marray,σarray)
end
function qnlab(j,s,nf,m,σ)
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   md = Int(2*m+(σtype(nf,σ)==2)+1)
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out[:,2] = abs.(out[:,2])
   out = kron(ones(Int,md),out)
   marray = kron(msbuilder(nf,m,σ),ones(Int,jsd))
   #[2j n ka kc m s]
   out = hcat(fill(Int(2*j),size(out,1)),out,marray,fill(σ,jsd*md))
   return out
end
#function qnlabv(j,s,nf,vtm,σ)
#   σt = σtype(nf,σ)
#   nlist = Δlist(j,s)
#   jsd = Int((2*j+1)*(2*s+1))
#   vd = Int(vtm+1)
#   out = zeros(Int,0,3)
#   for n in nlist
#      nd = Int(2*n+1)
#      part = zeros(Int,nd,3)
#      part[:,1] = fill(n,nd)
#      part[:,2] = collect(Int,-n:n)
#      part[:,3] = k2kc.(part[:,1],part[:,2])
#      out = vcat(out,part)
#   end
#   out[:,2] = abs.(out[:,2])
#   out = kron(ones(Int,vd),out)
#   vtrray = kron(nf .* vtcoll(vtm,σt) .+ σ,ones(Int,jsd))
#   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
#   return out
#end
function qnlab(j,s)
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      out = vcat(out,part)
   end
   out = kron(ones(Int,md),out)
   #[2j n ka kc m s]
   #out = hcat(fill(Int(2*j),size(out,1)),out,marray,fill(σ,jsd*md))
   return out
end
#function k2kc(n,k)
#"""
#Determines the value of Kc based on the value of N and |Kₐ|
#"""
#   ka = abs(k)
#   if k < 0
#      kc = n - ka + 1 - isodd(n + k)
#   elseif k == zero(k)
#      kc = n
#   else
#      kc = n - ka + isodd(n + k)
#   end
#   return kc
#end


function kgen(n::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,2*n+1)' )
end
function kgen2(n::Int)::LowerTriangular{Int, Matrix{Int}}
   ks = collect(-n:n)
   out = LowerTriangular(zeros(Int,length(ks),length(ks)))
   @simd for i in 1:length(ks)
      @inbounds out[i:end,i] = ks[i:end]
   end
   return out
end

function kgeni(n::Int,lb::Int)::Array{Int,2}
   return kron(collect(-n:n),ones(Int,lb)' )
end
function kgeni(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = collect(-n:n)
      si = fi + 1
   end
   return kron(out,ones(Int,lb)' )
end
function kgen(j::Float64,s::Float64)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,convert(Int,(2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = collect(-n:n)
      si = fi + 1
   end
   #lb = convert(Int,(2.0*j+1.0)*(2.0*s+1.0))
   return kron(out,ones(Int,convert(Int,(2.0*j+1.0)*(2.0*s+1.0)))' )
end
function ngen(n::Int)::Array{Int,2}
   return kron(fill(n,2*n+1),ones(Int,2*n+1)' )
end
function ngeni(j::Float64,s::Float64,lb::Int)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,Int((2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = fill(n,2*n+1)
      si = fi + 1
   end
   return kron(out,ones(Int,lb)' )
end
function ngen(j::Float64,s::Float64)::Array{Int,2}
   ns = Δlist(j,s)
   out = zeros(Int,convert(Int, (2.0*s+1.0)*(2.0*j+1.0)))
   si = 1
   for n in ns 
      fi = si + 2*n
      out[si:fi] = fill(n,2*n+1)
      si = fi + 1
   end
   return kron(out,ones(Int,convert(Int, (2.0*s+1.0)*(2.0*j+1.0)))' )
end

σcount(nfold::Int)::Int = floor(Int,nfold/2)+1
#function σtype(nfold,σ)
#   if σ==zero(σ) # A state
#      return 0
#   elseif (iseven(nfold))&&(σ==(σcount(nfold)-1)) # B state
#      return 2
#   else # E state
#      return 1
#   end
#end
#function msbuilder(T::Type,nfold::Number,mcalc::Number,σ::Number)
#   if nfold==0
#      return [1.]
#   else
#   σt = σtype(nfold,σ)
#   lim = mcalc*nfold
#   if σt==0
#      marray = collect(T,-lim:nfold:lim)
#   elseif σt==2
#      lim += σ
#      marray = collect(T,-lim:nfold:lim)
#   else
#      marray = collect(T,(-lim+σ):nfold:(lim+σ))
#   end
#   return marray
#   end
#end
msgen(mc::Int,nf::Int,σ::Int)::Array{Int} = msbuilder(Int,nf,mc,σ)
function mslimit(nfold,mcalc,σ)::Tuple{Int, Int}
   σt = σtype(nfold,σ)
   lim = mcalc*nfold
   if σt==0
      return -lim, lim
   elseif σt==2
      lim += σ
      return -lim, lim
   else
      return (-lim+σ), (lim+σ)
   end
end



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
   out[1] = (inp[3]+inp[1])/√2
   out[2] = inp[2]
   out[3] = -(inp[3]-inp[1])/√2
   return out
end
