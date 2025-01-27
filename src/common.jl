
################################################################################
#QN Related Functions
################################################################################

"""
This generates a list of all J₃ values that satisfy the Triangle Conditions with
   inputs J₁ and J₂. This is typically used to determine N values for a given
   J and S.
"""
function Δlist(J,S)
   return collect(Int(abs(J-S)):Int(J+S))
end

"""
This confirms that inputs J₁, J₂, and J₃ satisfy the Triangle Conditions with one another
"""
function Δtest(a,b,c)
   return c ⊆ collect(abs(a-b):abs(a+b))
end

"""
This function determines the signed K values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group for a given N.
   The arrays are returned as a tuple of A₁ then A₂
"""
function klist(n)
   Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
   Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
   if iseven(n)
      return Γ1, Γ2
   else
      return Γ2, Γ1
   end
end
"""
This function determines the signed Kₐ values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group for a given N and m pair.
   The arrays are returned as a tuple of A₁ then A₂. This should only be used
   for the A states of a G₆ molecule.
"""
function klist(n,m)
   Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
   Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
   if iseven(n)&&(m≥0)||isodd(n)&&(m<0)
      return Γ1, Γ2
   else
      return Γ2, Γ1
   end
end
"""
This function determines the signed Kₐ values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group across an entire -m:m range
   a given N.
   The arrays are returned as a tuple of A₁ then A₂.
"""
function mklist(n,mcalc)
   A1, A2 = klist(n,-mcalc)
   nd = convert(Int,2*n+1)
   shift = nd
   for m ∈ (1 - mcalc):mcalc
      a1, a2 = klist(n,m)
      A1 = vcat(A1,a1 .+ shift)
      A2 = vcat(A2,a2 .+ shift)
      shift += nd
   end
   return A1, A2
end

"""
Determines the value of Kc based on the value of N and |Kₐ|
"""
function k2kc(n,k)
   ka = abs(k)
   if k < 0
      kc = n - ka + 1 - isodd(n + k)
   elseif k == zero(k)
      kc = n
   else
      kc = n - ka + isodd(n + k)
   end
   return kc
end
function kakc2k(n,ka,kc)
   ka = abs(ka)
   if iseven(n+ka+kc)
      return ka*powneg1(n+ka)
   else isodd(n+ka+kc)
      return -ka*powneg1(n+ka)
   end
end

"""
This determines the specific index of a N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
function qn2ind(n,ka,kc)
   ka = kakc2k(n,ka,kc)
   ind = convert(Int,n*(n+1) + ka +1)
end
"""
This determines the specific index of a J, S, N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
function qn2ind(j,s,n,ka,kc)
   jp = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)
   np = sum(2 .* collect((j-s):(n-1)) .+ 1)
   kp = n + kakc2k(n,ka,kc) + 1
   ind = jp + np + kp
   ind = convert(Int,ind)
   return ind
end
"""
This determines the specific index of a J, S, N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
   This is the old version based on the torsional
"""
#function qn2indm(nf,mcalc,m,j,s,n,ka,kc)
#   if nf≠zero(nf)
#   #ka = abs(ka)
#   ind = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)*(2*mcalc+1)
#   ind += (mcalc+floor(m/nf))*(2*s+1)*(2*j+1)
#   ind += sum(2 .* collect((j-s):(n-1)) .+ 1) + n + kakc2k(n,ka,kc) + 1
#   ind = convert(Int,ind)
#   return ind
#   else
#   return qn2ind(j,s,n,ka,kc)
#   end
#end
"""
This determines the specific index of a J, S, N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
function qn2ind_prev(nf,vtm,m,j,s,n,ka,kc)
   if (nf≠zero(nf))&&(vtm ≠ 0)
   #ka = abs(ka)
   ind = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)*(2*mcalc+1)
   ind += (mcalc+floor(m/nf))*(2*s+1)*(2*j+1)
   ind += sum(2 .* collect((j-s):(n-1)) .+ 1) + n + kakc2k(n,ka,kc) + 1
   ind = convert(Int,ind)
   return ind
   else
   return qn2ind(j,s,n,ka,kc)
   end
end
function qn2ind(nf,vtm,m,j,s,n,ka,kc)
   if (nf≠zero(nf))
   #ka = abs(ka)
   ind = sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)*(vtm+1)
   #ind += (vtm+floor(Int,m/nf))*(2*j+1) #<--- This function is wrong for vt≠0
   ind *= Int(2s+1)
   ind += (ceil(Int,0.5*vtm) + floor(Int,m/nf))*(2*j+1)
   ind += sum(2 .* collect((j-s):(n-1)) .+ 1) + n + kakc2k(n,ka,kc) + 1
   ind = convert(Int,ind)
   return ind
   else
   return qn2ind(j,s,n,ka,kc)
   end
end

"""
Determines the number of σ values that will result in separable torsional blocks a provided
   nfold value. So for nfold = 3, it returns 2 since we have the A states of
   σ=0 and E states of σ=1. It will also return 2 for nfold = 2 since we have A
   and B states.
"""
function σcount(nfold::Real)::Int
   if isodd(nfold)
      out = ceil(Int,0.5*nfold)
   elseif iseven(nfold)&&nfold≠0
      out = floor(Int,nfold/4) + 1
   else #nfold == 0
      out = 1
   end
   return out
end

function msgen(T::Type,nfold::Real,mcalc::Real,σ::Real)
   if nfold==0
      return ones(T,1)
   else
   lim = mcalc*nfold
   if (σ==0)&&(isodd(nfold))
      marray = collect(T,-lim:nfold:lim)
   elseif (σ==0)&&(iseven(nfold))
      lim = floor(Int,lim/2)
      marray = collect(T,-lim:floor(T,nfold/2):lim)
   elseif (σ≠0)&&(iseven(nfold))
      lim = floor(Int,lim/2)
      marray = collect(T,-lim+σ:floor(T,nfold/2):lim+σ)
   else
      marray = collect(T,(-lim+σ):nfold:(lim+σ))
   end
   return marray
   end
end
function msgen(nfold::Int,mcalc::Int,σ::Int)
   marray = msgen(Int,nfold,mcalc,abs(σ))
   if σ < 0
      marray .*= -1
   end
   return marray
end

"""
This calculates a lot of values helpful for the spin-rotation matrix structure:
   the list of n values, the number of K states for each n, an array of the
   start and end indices for each n in the list, and the number of substates for
   the given J and S pair. The array of indices has starts on the [:,1] values
   and the ends on the [:,2] values.
"""
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
"""
This expands the srprep function to the spin-torsion-rotation matrix structure
   by including the number of m values, md = 2*mcalc+1
"""
function srprep(J,S,md)
   ns = Δlist(J,S)
   nd = 2 .* Int.(ns) .+ 1
   out = ones(Int, length(ns),2)
   out[1,2] = nd[1]*md
   for i in 2:length(ns)
      out[i,1] = out[i-1,2] + 1
      out[i,2] = out[i,1] + md*nd[i] - 1
   end
   jd = Int((2*S+1)*(2*J+1))
   return ns, nd, out, jd
end

"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
function jinds(j,s)
   snd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)) +1
   fnd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):j) .+ 1))
   return snd,fnd
end
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
function jlinds(j,s)
   snd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)) +1
   fnd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):j) .+ 1))
   return collect(snd:fnd)
end

function jinds(j,s,m)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1)*(2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)) +1
   fnd = convert(Int, (2*m+1)*(2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):j) .+ 1))
   return snd,fnd
end
function jlinds(j,s,m)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*m+1)*(2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)) +1
   fnd = convert(Int, (2*m+1)*(2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):j) .+ 1))
   return collect(snd:fnd)
end
function mfilter(jsd,nf,mc,mm,σ)
   mlowest = σ
   mhighst = (σ + nf*mm)*powneg1(isodd(mm))
end

function cart2sphr(inp::Array{Float64,1})::Array{Float64,1}
   #assumes a [z x y] ordering
   out = zeros(3)
   out[1] = (inp[2]+inp[3])/√2
   out[2] = inp[1]
   out[3] = -(inp[2]-inp[3])/√2
   return out
end

################################################################################
#These are the various Wang Matrices
################################################################################
function eye(x::Int)::SparseMatrixCSC{Float64, Int64}
"""
Returns a sparse identity matrix of size x. Mostly used for type consistency in
   the Wang Transformations.
"""
   return sparse(I,x,x)
end
function un(m::Int)::SparseMatrixCSC{Float64, Int64}
"""
This generates a standard Wang Transformation Matrix for a matrix of size m. As
   an example, if you want the matrix for N=1, input m = 3 (= 2*N+1). Or if you
   want to span the even values of the free rotor from -4:4, such as for a 2-fold
   rotor's B states, input m = 4. Even m values will exclude the untransformed
   0 state that centers the typical form of these matrices.
"""
   if iseven(Int(m))
      n = Int(m/2)
      out = (1/sqrt(2)) .* [-eye(n) rotl90(eye(n)); rotl90(eye(n)) eye(n)]
      return out
   elseif isodd(Int(m))
      n = Int((m-1)/2)
      out = (sqrt(0.5)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
       rotl90(eye(n)) zeros(n) eye(n)]
      return out
   end
end

function ut(m,σt)::SparseMatrixCSC{Float64, Int64}
"""
This builds the torsional Wang Transformation matrix for a span of -m:m with a
   rotational block size of 2*n+1. A pure torsional form can be built using n=0
"""
   if σt == 0
      out = Diagonal(vcat(fill(-√.5,m), 1.0, fill(√.5,m)))
      out += rotl90(Diagonal(vcat(fill(√.5,m), 0.0, fill(√.5,m))))
      out = sparse(out)
   elseif σt==2
      m += 1
      out = Diagonal(vcat(fill(-√.5,m), fill(√.5,m)))
      out += rotl90(Diagonal(vcat(fill(√.5,m), fill(√.5,m))))
      out = sparse(out)
   else
      out = sparse(I,2*m+1,2*m+1)
   end
   return out
end

function eyr(x::Int)::Array{Float64,2}
   diagm(ones(x))
end
function ur(n::Int)::SparseMatrixCSC{Float64, Int64}
   out = spdiagm(vcat(fill(-√.5,n), 1.0, fill(√.5,n)))
   out += rotl90(spdiagm(vcat(fill(√.5,n), 0.0, fill(√.5,n))))
   return out
end
function ur(j::Float64,s::Float64)::SparseMatrixCSC{Float64, Int64}
"""
This builds the rotational Wang Transformation matrix for every n in Δlist(j,s).
   This will be kronecker producted with an identy matrix of size 2*m+1 for the
   torsional-rotational nature. A purely rotational form can be built using m=0
"""
   ns, nd, ni, jsd = srprep(j,s)
   out = spzeros(jsd,jsd)
   for i in 1:length(ns)
      r = ni[i,1]:ni[i,2]
      out[r, r] = ur(ns[i])
   end
   return out
end

vt2m(vtm)::Int = ceil(vtm/2)*cospi(vtm)
#vt2m(vtm,σt)::Int = ceil(vtm/2)*cospi(vtm) + (vtm≤zero(vtm))*δi(σt,2)
function vtlist(vtm) 
   if vtm == zero(vtm)
      return [vt2m(vtm)]
   else
      return sort([vt2m(vtm); vt2m(vtm-1)])
   end
end
function vtcoll(vtm) 
   if vtm == zero(vtm)
      return [vt2m(vtm)]
   else
      a = vt2m(vtm)
      b = vt2m(vtm-1)
      ref = sort([a; b])
      return collect(ref[1]:ref[2])
   end
end


################################################################################
# Parameter shifting functions
################################################################################

function PAM2RAM(inps)
"""
Hypothetically this rotates the parameters from the Principal Axis values
   to the Rho Axis values. Does it work? Maybe! I'm not super sure at this moment
interal parameters:
     1  2  3   4   5   6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; Fρ; V3; ao; a; b; d; η]
input parameters:
     1  2  3  4  5   6   7   8   9   10  11
prs =[A; B; C; δ; F; V3; ϵzz; ϵxx; ϵyy; ϵxzb; η]
"""
   out = zeros(Float64,12)
   cδ = cos(inps[4])
   sδ = sin(inps[4])
   Ap = inps[1]
   Bp = inps[2]
   F0 = inps[5]
   Ar = Ap*Bp/(Bp*cδ^2 + Ap*sδ^2)
   Br = Ap*Bp/(Ap*cδ^2 + Bp*sδ^2)
   Fr = F0^2/(F0-Ar)
   Aeff = Ar + Ar^2/(F0-Ar)
   out[1] = Aeff
   out[2] = Br
   out[3] = inps[3]  #Cr = Cp
   out[4] = cδ*sδ*(1/Bp - 1/Ap)/(2*Ar*Br)
   out[5] = Fr
   out[6] = F0*Ar/(F0-Ar)
   out[7] = inps[6] #V3
   out[8] = -(inps[7]+inps[9]+inps[8])/3.0
   out[9] = -(2.0*inps[7]-inps[9]-inps[8])/6.0
   out[10] = (inps[8] - inps[9])*0.5
   out[11] = inps[10] #ϵxzb
   out[12] = inps[11] #η
   return out
end
function RAM2PAM(ints)
"""
Hypothetically this rotates the parameters from the Rho Axis values to the
   Principal Axis values. Does it work? Maybe! I'm not super sure at this moment
interal parameters:
     1  2  3   4   5   6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; Fρ; V3; ao; a; b; d; η]
input parameters:
     1  2  3  4  5   6   7   8   9   10  11
prs =[A; B; C; δ; F; V3; ϵzz; ϵxx; ϵyy; ϵxzb; η]
"""
   out = zeros(Float64,11)
   Fr = ints[5]
   ρ = ints[6]/ints[5]
   Ar = ints[1] - Fr*ρ^2
   Br = ints[2]
   Dab = ints[4]
   F0 = Ar/ρ
   Bp, Ap = eigvals([Ar Dab; Dab Br])
   out[1] = Ap
   out[2] = Bp
   out[3] = ints[3]
   out[4] = 0.5*asin(4*Ar*Br*Dab/(1/Bp - 1/Ap))
   out[5] = F0
   out[6] = ints[7]
   a0 = ints[8]
   ϵzz = -(a0 + 2*ints[9])
   out[7] = ϵzz
   out[8] = b - 0.5*(ϵzz + 3*a0)
   out[9] = -b - 0.5*(ϵzz + 3*a0)
   out[10] = ints[11]
   out[11] = ints[12]
   return out
end
function paramprep(params)
"""
As I am growingly convinced that ρ is best treated as a derived parameter instead
   of a independent one, this determines ρ from A and F, reduces F into Fr,
   calculates ρFr, and changes A_RAM into A_eff = A_RAM + Fr*ρ^2. The shift is
   done for the ease of the analytic derivatives
input parameters:
      1  2  3   4   5   6   7  8  9 10 11 12
prs =[A; B; C; Dab; F; V3; ao; a; b; d; η]
internal parameters:
      1  2  3   4    5   6   7   8  9 10 11 12
prs =[A; B; C; Dab; Fr; Fρ; V3; ao; a; b; d; η]
"""
   A = params[1]
   F = params[5]
   Fr = F^2/(F-A)
   ρ = A/F
   Fρ = Fr*ρ
   ashift = Fρ*ρ
   output = zeros(13)
   output[1] = A + ashift
   output[2] = params[2]
   output[3] = params[3]
   output[4] = params[4]
   output[5] = Fr
   output[6] = Fρ
   output[7] = params[6]
   output[8] = params[7]
   output[9] = params[8]
   output[10] = params[9]
   output[11] = params[10]
   output[12] = params[11]
   output[13] = params[12]
   return output
end

function greaterof(x::Real,y::Real)
   if x ≥ y
      return x
   else
      return y
   end
end
function lesserof(x::Real,y::Real)
   if x ≥ y
      return y
   else
      return x
   end
end

function sfdiv(a::Real,b::Real)::Float64
   if b ≤ zero(b)
      return 0.0
   else
      return a / b
   end
end