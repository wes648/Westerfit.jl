function Rotbuilder(pr,nmax)
   nmd = convert(Int,2*nmax+1)
   out = zeros(typeof(pr[1]),nmd,nmd,nmax+1)
   @threads for n in 1:nmax
      nd = 2*n+1
      @inbounds out[nd,nd,n+1] = Hrot(pr,n)
   end
   return out
end
function Rotcall(n::Int,bigrot)
   nd = 2*n+1
   return bigrot[nd,nd,n+1]
end
function Spibuilder(pr,hrot,jlist,s)
   sd = convert(Int,2.0*s+1.0)
   jd = convert(Int,2.0*jlist[end]+1.0)*sd
   out = spzeros(jd,jd,length(jlist))
   @threads for i in 1:length(jlist)
      j = jlist[i]
      jd = convert(Int,2.0*j+1.0)*sd
      @inbounds out[jd,jd,i] = Hsr2(pr,hrot,j,s)
      @inbounds out[jd,jd,i] = Hhyp!(pr,out[jd,jd,i],j,s)
   end
   return out
end
function Spicall(j,s,bigsr)
   ji = convert(Int,j-0.5*iseven(Int(2.0*s+1.0))) + 1
   jd = convert(Int,(2.0*j+1.0)*(2.0*s+1.0))
   return bigsr[jd,jd,ji]
end
function torbuild(pr,hsr,mcalc,j,s,σ)
   out = kron(I(Int(2*mcalc+1)), Spicall(j,s,hsr))
   out += Htor(pr,mcalc,j,s,σ)
   return out
end
function tsrdiag_build(pr,j,s,refmat,mcalc,σ)
"""
Builds, diagonalizes, and assigns the Hamiltonian for a given J, S, σ, & mcalc
"""
   H = Spicall(j,s,refmat)
   U = ur(j,s,mcalc)
   if σ==zero(σ)
      U *= ut(mcalc,j,s)
      H = Matrix(U*H*U)
      #H = Matrix(Htsrmat2(pr,j,s,mcalc,σ))
   else
      H = Matrix(U*H*U)
   end
   if j≥1.0
      H, rvec = jacobisweep2(H,Int(floor((j+s+6*mcalc)/8)))
      vals, vecs = LAPACK.syev!('V', 'U', H)
   else
      vals, vecs = LAPACK.syev!('V', 'U', H)
      rvec = I(length(vals))
   end
   qns, vals, vecs = assign(j,s,σ,mcalc,vals,vecs,rvec)
   return qns, vals, vecs
end

function tsrcalc_build(prm,s,nmax,mcalc,mmax,σ)
"""
Calculates all the energy levels, eigenvectors, & quantum numbers for all J values
   for a given σ
"""
   refmat = Rotbuilder(prm,nmax)
   println("Rot mats built!")
   jmin =  0.5*iseven(Int(2*s+1))
   jmax = nmax - s
   js = collect(Float64,jmin:jmax)
   refmat = Spibuilder(prm,refmat,js,s)
   println("Spi mats built!")
   σs = σcount(NFOLD)
   jsσs = zeros(0,2)
   for σ in 1:σ2
      jsσs = vcat(jsσs, [jarray fill(σ-1,length(js))])
   end
   jd = Int((2*s+1)*sum(2.0 .* js .+ 1.0)) #This looks WAY too big -W 9/20/22
   outvals = zeros(Float64,Int(jd*(2*mmax+1)),σs)
   outvecs = zeros(Float64,Int((2*s+1)*(2*jmax+2)*(2*mcalc+1)),Int(jd*(2*mmax+1)),σs)
   outqns = zeros(Float64,Int(jd*(2*mmax+1)),6,σs)
   @threads for i in 1:size(jsσs)[1]
      j = jsσs[i,1]
      σ = jsσs[i,2]
      sind, find = jinds(j,s,mmax)
      tqns, tvals, tvecs = tsrdiag_build(prm,j,s,refmat,mcalc,σ)
      si = findfirst(isequal(-mmax*NFOLD+σ),tqns[:,5])
      fi = findlast(isequal(mmax*NFOLD+σ),tqns[:,5])
      outvals[sind:find,σ+1] = tvals[si:fi]
      outvecs[1:Int((2*s+1)*(2*j+1)*(2*mcalc+1)),sind:find,σ+1] = tvecs[:,si:fi]
      outqns[sind:find,:,σ+1] = tqns[si:fi,:]
   end
   return outqns, outvals, outvecs
end
