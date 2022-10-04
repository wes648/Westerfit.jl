
################################################################################
#QN Related Functions
################################################################################

function Δlist(J,S)
"""
This generates a list of all J₃ values that satisfy the Triangle Conditions with
   inputs J₁ and J₂. This is typically used to determine N values for a given
   J and S.
"""
   max = Int(J+S)
   min = Int(abs(J-S))
   return collect(min:max)
end

function Δtest(a,b,c)
"""
This confirms that inputs J₁, J₂, and J₃ satisfy the Triangle Conditions with one another
"""
   return c ⊆ collect(abs(a-b):abs(a+b))
end

function klist(n)
"""
This function determines the signed K values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group for a given N.
   The arrays are returned as a tuple of A₁ then A₂
"""
   Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
   Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
   if iseven(n)
      return Γ1, Γ2
   else
      return Γ2, Γ1
   end
end
function klist(n,m)
"""
This function determines the signed Kₐ values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group for a given N and m pair.
   The arrays are returned as a tuple of A₁ then A₂. This should only be used
   for the A states of a G₆ molecule.
"""
   Γ1 = sort([collect(-1:-2:-n); collect(0:2:n)]) .+ (n + 1)
   Γ2 = sort([collect(-2:-2:-n); collect(1:2:n)]) .+ (n + 1)
   if iseven(n)&&(m≥0)||isodd(n)&&(m<0)
      return Γ1, Γ2
   else
      return Γ2, Γ1
   end
end
function mklist(n,mcalc)
"""
This function determines the signed Kₐ values that correspond to the A₁ and A₂
   irreducible representations of the Cₛ point group across an entire -m:m range
   a given N.
   The arrays are returned as a tuple of A₁ then A₂.
"""
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

function qngen(n,m,σ)
"""
This generates a list of all the torsion-rotation quantum numbers across a
   -m:m range for a given N and σ pair. Be aware that this uses the signed Kₐ
   convention. Kc is provided for user convenience but not generally used by
   the program as the sign of Kₐ accounts for this.The output columns are
   ordered as N, Kₐ, Kc, m, σ
"""
   nd = Int(2*n+1)
   md = Int(2*m+1)
   narray = fill(n,nd*md)
   karray = kron(ones(Int,md),collect(Int,-n:n))
   kcarray = k2kc.(narray,karray)
   marray = kron(NFOLD .* collect(Int,-m:m) .+ σ,ones(Int,nd))
   σarray = fill(σ,nd*md)
   out = hcat(narray,karray,kcarray,marray,σarray)
end
function qngen(j,s,m,σ)
"""
This generates a list of all the spin-torsion-rotation quantum numbers across a
   -m:m range for a given J, S, and σ set. Be aware that this uses the signed Kₐ
   convention. Kc is provided for user convenience but not generally used by
   the program as the sign of Kₐ accounts for this. The output columns are
   ordered as J, N, Kₐ, Kc, m, σ. S is not included under the hope that it is
   recorded by the user elsewhere
"""
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   md = Int(2*m+1)
   out = zeros(0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out = kron(ones(Int,md),out)
   ms = kron(NFOLD .* collect(Int,-m:m) .+ σ,ones(Int,jsd))
   out = hcat(fill(j,size(out)[1]),out,ms,fill(σ,jsd*md))
   return out
end

function k2kc(n,k)
"""
Determines the value of Kc based on the value of N and |Kₐ|
"""
   ka = abs(k)
   if k < 0
      kc = n - ka + 1
   else
      kc = n - ka
   end
   return kc
end

function qn2ind(n,ka,kc)
"""
This determines the specific index of a N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
   ka = abs(ka)
   ind = convert(Int,n*(n+1)+ka*(-1)^(n-ka-kc) +1)
end
function qn2ind(j,s,n,ka,kc)
"""
This determines the specific index of a J, S, N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
   ka = abs(ka)
   jp = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)
   np = sum(2 .* collect((j-s):(n-1)) .+ 1)
   kp = n + ka*(-1)^(n-ka-kc) + 1
   ind = jp + np + kp
   ind = convert(Int,ind)
   return ind
end
function qn2ind(mcalc,m,j,s,n,ka,kc)
"""
This determines the specific index of a J, S, N, Ka, Kc state in the large array.
   Ka is assumed to not be signed for better compatibility with literature.
"""
   if NFOLD!=zero(NFOLD)
   ka = abs(ka)
   ind = (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)*(2*mcalc+1)
   ind += (mcalc+floor(m/NFOLD))*(2*s+1)*(2*j+1)
   ind += sum(2 .* collect((j-s):(n-1)) .+ 1) + n + ka*(-1)^(n-ka-kc) + 1
   ind = convert(Int,ind)
   return ind
   else
   return qn2ind(j,s,n,ka,kc)
   end
end

function σcount(nfold)
"""
Determines the number of σ values that will result in unique values a provided
   nfold value. So for nfold = 3, it returns 2 since we have the A states of
   σ=0 and E states of σ=1. It will also return 2 for nfold = 2 since we have A
   and B states.
"""
   out = floor(Int,nfold/2)+1
end

function msetter(nfold,mcalc,mmax)
"""
Sets mcalc & mmax to zero if NFOLD=0 to fully disable torsional behavior.
   If NFOLD>0, makes sure that mcalc is not less than NFOLD. This is due to the
   odd behavior of westerfit which will fail if the torsional basis is not
   sufficiently large. While technically a bug, there are no current plans to
   change this behavior.
"""
   if nfold==zero(nfold)
      mcalc = zero(mcalc)
      mmax = zero(mmax)
   else
      if mcalc < nfold
         println("mcalc was too low for NFOLD value. Setting mcalc=NFOLD")
         mcalc = nfold
      end
   end
   return mcalc, mmax
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
function srprep(J,S,md)
"""
This expands the srprep function to the spin-torsion-rotation matrix structure
   by including the number of m values, md = 2*mcalc+1
"""
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

function jinds(j,s)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):(j-1)) .+ 1)) +1
   fnd = convert(Int, (2*s+1)*sum(2 .* collect((0.5*isodd(2*s)):j) .+ 1))
   return snd,fnd
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

################################################################################
#These are the various Wang Matrices
################################################################################
function eye(x)
"""
Returns a sparse identity matrix of size x. Mostly used for type consistency in
   the Wang Transformations.
"""
   spdiagm(ones(x))
end
function un(m)
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
      out = (1/sqrt(2)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
       rotl90(eye(n)) zeros(n) eye(n)]
      return out
   end
end

function ut(m,n)
"""
This builds the torsional Wang Transformation matrix for a span of -m:m with a
   rotational block size of 2*n+1. A pure torsional form can be built using n=0
"""
   nd = 2*n+1
   out = (1/sqrt(2)) .*  [-1*eye(m) zeros(m) rotl90(eye(m)); zeros(1,m) sqrt(2) zeros(1,m);
   rotl90(eye(m)) zeros(m) eye(m)]
   out = kron(out,eye(nd))
   return out
end
function ut(m,j,s)
"""
This builds the torsional Wang Transformation matrix for a span of -m:m with
   rotational blocks for every n in Δlist(j,s)
"""
   jd = Int((2*s+1)*(2*j+1))
   out = (1/sqrt(2)) .*  [-1*eye(m) zeros(m) rotl90(eye(m)); zeros(1,m) sqrt(2) zeros(1,m);
   rotl90(eye(m)) zeros(m) eye(m)]
   out = kron(out,eye(jd))
   return out
end
function ur(n,m)
"""
This builds the rotational Wang Transformation matrix for a given n. This will
   be kronecker producted with an identy matrix of size 2*m+1 for the
   torsional-rotational nature. A purely rotational form can be built using m=0
"""
   md = 2*m+1
   out = (1/sqrt(2)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
      rotl90(eye(n)) zeros(n) eye(n)]
   out = kron(eye(md),out)
   return out
end
function ur(j,s,m)
"""
This builds the rotational Wang Transformation matrix for every n in Δlist(j,s).
   This will be kronecker producted with an identy matrix of size 2*m+1 for the
   torsional-rotational nature. A purely rotational form can be built using m=0
"""
   nlist = Δlist(j,s)
   out = zeros(0,0)
   for n in nlist
      part = (1/sqrt(2)) .* [-eye(n) zeros(n) rotl90(eye(n)); zeros(1,n) sqrt(2) zeros(1,n);
         rotl90(eye(n)) zeros(n) eye(n)]
      out = cat(out,part,dims=(1,2))
   end
   md = 2*m+1
   out = kron(eye(md),out)
   return out
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

function greaterof(x::Number,y::Number)
   if x ≥ y
      return x
   else
      return y
   end
end
function lesserof(x::Number,y::Number)
   if x ≥ y
      return y
   else
      return x
   end
end
