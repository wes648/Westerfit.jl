"""
Jacobi transformers to for cleaner eigenvectors.
Blatantly stolen & modified from
https://github.com/albi3ro/M4/blob/master/Numerics_Prog/Jacobi-Transformation.ipynb
Thanks Christina Lee!
I have removed the compV options since I always want to do that in the westerfit package
"""

#module jacobi
#
#export jacobisweep

using LinearAlgebra, LinearAlgebra.BLAS

# First, Lets make our nice, helpful functions

## A function to look at the convergence
function convergence(A::Array)
   num = 0.0
   l=size(A,1)
   for ii in 1:(l-1)
      for jj in (ii+1):l ## just looking at the lower triangle
         num += A[ii,jj]^2
         #println(ii,' ',jj,' ',num,' ',A[ii,jj])
      end
   end
   return num
end

# This makes a matrix easier to look at when its filled
# with 1.043848974e-12 everywhere
function roundmatrix(A::Array,rtol::Real)
   Ap=copy(A)
   for ii in 1:length(A)
      if abs(Ap[ii]) < rtol
         Ap[ii] = 0.0
      end
   end
   return Ap;
end

function makeA(n)
   A=randn(n,n);
   for ii in 1:n
      A[ii,1:ii]=transpose(A[1:ii,ii])
   end
   V=Matrix{Float64}(I,n,n) #initializing the orthogonal transformation
   return A#,copy(A),V
end

function nulldiag(A::Array)::Array
   mat = copy(A)
   mat[diagind(A)] .= 0
   return mat
end
function lrgstoffd(A::Array)
   mat = nulldiag(A)
   out = [iamax(mat[:,i]) for i in 1:size(mat,2)]
   return out
end

#function eye(x::Int)::SparseMatrixCSC{Float64, Int64}
#"""
#Returns a sparse identity matrix of size x. Mostly used for type consistency in
#   the Wang Transformations.
#"""
#   return sparse!(I,x,x)
#end
function givenS(A,p::Int,q::Int,θ::Float64)
   V = eye(size(A,1))
   V[p,p] = cos(θ)
   V[q,q] = cos(θ)
   V[p,q] = sin(θ)
   V[q,p] = -sin(θ)
   return V
end
#Now on to the Rotations!

# We don't always want to compute the eigenvectors, so those are in the
# optional entries slot.
# Both tell the function to compute the vectors with computeV=true
# and input the V=V after the semicolon.
function RotateFaster(A,p::Int,q::Int)
   Apq = A[p,q]
   Δpq = A[q,q] - A[p,p]
   if (Apq != 0.0)&&(Δpq != 0.0)
      θ = 0.5*atan(2.0*Apq, abs(Δpq))*sign(Δpq)
      G = LinearAlgebra.Givens(p,q,cos(θ),sin(θ))
      A = dropzeros!(G' * A * G)
   else
      G = LinearAlgebra.Givens(p,q,1.0,0.0)
   end
   return A, G
end

function Rotate(A::Array,p::Int,q::Int, V::Array=Matrix{Float64}(I,1,1) )
   Apq = A[p,q]
   if Apq != 0.0
   θ = (A[q,q] - A[p,p]) / (2.0*A[p,q])
   t = sign(θ)/(abs(θ)+sqrt(θ^2 + 1.0))
   c = 1.0/√(t^2 + 1.0)
   s = t*c
   τ = s/(1.0 + c)
   l=size(A,1)
   Ap=copy(A[:,p])
   Aq=copy(A[:,q])
   rs = collect(1:l)
   @. A[rs,p] = Ap[rs] - s*(Aq[rs] + τ*Ap[rs])
   @. A[rs,q] = Aq[rs] + s*(Ap[rs] - τ*Aq[rs])
   @. A[p,rs] = A[rs,p]
   @. A[q,rs] = A[rs,q]
   A[p,q] = 0.0
   A[q,p] = 0.0
   A[p,p] = Ap[p] - t*Aq[p]
   A[q,q] = Aq[q] + t*Aq[p]
   Vp = copy(V[:,p])
   Vq = copy(V[:,q])
   @. V[rs,p] = c*Vp[rs] - s*Vq[rs]
   @. V[rs,q] = s*Vp[rs] + c*Vq[rs]
   end
   return A,V
end
function Rotate(A::SparseMatrixCSC{Float64,Int64},p::Int,q::Int, V)
   if A[p,q] != 0.0
   t = (A[q,q] - A[p,p]) / (2.0*A[p,q]) #this is the angle, called t to save mem
   t = sign(t)/(abs(t)+sqrt(t^2 + 1.0))
   c = 1.0/√(t^2 + 1.0)
   s = t*c
   τ = s/(1.0 + c)
   Ap=copy(A[:,p])
   Aq=copy(A[:,q])
   rs = collect(1:size(A,1))
   @. A[rs,p] = Ap[rs] - s*(Aq[rs] + τ*Ap[rs])
   @. A[rs,q] = Aq[rs] + s*(Ap[rs] - τ*Aq[rs])
   @. A[p,rs]=A[rs,p]
   @. A[q,rs]=A[rs,q]
   A[p,q] = 0.0
   A[q,p] = 0.0
   A[p,p] = Ap[p] - t*Aq[p]
   A[q,q] = Aq[q] + t*Aq[p]
   Vp = copy(V[:,p])
   Vq = copy(V[:,q])
   @. V[rs,p] = c*Vp[rs] - s*Vq[rs]
   @. V[rs,q] = s*Vp[rs] + c*Vq[rs]
   end
   return A,V
end
# This function performs one sweep
function Sweep(A,V=Matrix{Float64}(I,size(A)))
   for ii in 2:size(A,1)
      for jj in 1:(ii-1) ## Just over one triangle
         A, V = Rotate(A,ii,jj,V)
      end
   end
   return A,V
end
function SweepF(A)
   V = eye(size(A,1))
   for ii in 2:size(A,1)
      for jj in 1:(ii-1) ## Just over one triangle
      A, G = RotateFaster(A,ii,jj)
      V = G*V
      end
   end
   return A,V
end
function SparseSweep(A,V=Matrix{Float64}(I,size(A)))
   r, c, = findnz(A)
   rf = r[c .< r]
   cf = c[c .< r]
   for i in 1:length(rf)
      A, V = Rotate(A,rf[i],cf[i],V)
   end
   return A,V
end
function SparseSweepF(A,V = eye(size(A,1)), f=0)
   r, c, = findnz(A)
   #perm = [c .< (r .- f)]
   rf = r[c .< (r .- f)]
   cf = c[c .< (r .- f)]
   for i in 1:length(rf)
      A, G = RotateFaster(A,rf[i],cf[i])
      V = G*V
   end
   return A,V
end
function jacobisparse!(A,cnt=1,f=0)
   V = eye(size(A,1))
   for i in 1:cnt 
      A,V = SparseSweepF(A,V,f)
   end
   return A, V
end

function SparseList(A)
   r, c, = findnz(A)
   rf = r[c .< r]
   cf = c[c .< r]
   for i in 1:length(rf)
      println([rf[i] cf[i]])
   end
end

function limsparsweep(A,cnt=1)
   V = eye(size(A,1))
   for m in 1:cnt
      A,V = SparseSweep(A,V)
   end
   return A,V
end

function limsweep(A,cnt=1,V=Matrix{Float64}(I,size(A)))
   n=size(A,1)
   for m in 1:cnt
   for ii in 2:n
      jjlist = lrgstoffd(A)
      for jj in jjlist 
         A, V = Rotate(A,ii,jj,V)
      end
   end
   end
   return A,V
end


function jacobisweep(A,iters)
   n=size(A,1)
   V = Matrix(Float64.(I(n)))
   Vout = Matrix(Float64.(I(n)))
   for i in 1:iters
      A,V = Sweep(A,V)
      Vout = transpose(V)*Vout
   end
   return Matrix(A),Vout
end
function jacobisweep2(A,iters)
   n=size(A,1)
   V = Matrix(Float64.(I(n)))
   Vout = Matrix(Float64.(I(n)))
   for i in 1:iters
      A,V = Sweep(A,V)
#      Vout = transpose(V)*Vout
   end
   return Matrix(A),V
end
#end
