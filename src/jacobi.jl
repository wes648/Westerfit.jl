"""
Jacobi transformers to for cleaner eigenvectors.
Blatantly stolen & modified from
https://github.com/albi3ro/M4/blob/master/Numerics_Prog/Jacobi-Transformation.ipynb
Thanks Christina Lee!
"""

#module jacobi
#
#export jacobisweep

using LinearAlgebra

# First, Lets make our nice, helpful functions

## A function to look at the convergence
function convergence(A::Array)
    num=0.0
    l=size(A)[1]
    for ii in 1:(l-1)
        for jj in (ii+1):l ## just looking at the lower triangle
            num+=A[ii,jj]^2
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
        if abs(Ap[ii])<rtol
            Ap[ii]=0
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
    return A,copy(A),V
end

#Now on to the Rotations!

# We don't always want to compute the eigenvectors, so those are in the
# optional entries slot.
# Both tell the function to compute the vectors with computeV=true
# and input the V=V after the semicolon.

function Rotate(A::Array,p::Int,q::Int; computeV=false, V::Array=Matrix{Float64}(I,1,1) )
    θ=(A[q,q]-A[p,p])/(2*A[p,q]);
    t=sign(θ)/(abs(θ)+sqrt(θ^2+1));
    c=1/sqrt(t^2+1)
    s=t*c
    τ=s/(1+c)
    l=size(A)[1]
    Ap=copy(A[:,p])
    Aq=copy(A[:,q])
    for r in 1:l
        A[r,p]=Ap[r]-s*(Aq[r]+τ*Ap[r])
        A[r,q]=Aq[r]+s*(Ap[r]-τ*Aq[r])
        A[p,r]=A[r,p]
        A[q,r]=A[r,q]
    end
    A[p,q]=0
    A[q,p]=0
    A[p,p]=Ap[p]-t*Aq[p]
    A[q,q]=Aq[q]+t*Aq[p]
    if computeV==true
        Vp=copy(V[:,p])
        Vq=copy(V[:,q])
        for r in 1:l
            V[r,p]=c*Vp[r]-s*Vq[r]
            V[r,q]=s*Vp[r]+c*Vq[r]
        end
        return A,V
    else
        return A;
    end
end

# This function performs one sweep
function Sweep(A;compV=false,V=Matrix{Float64}(I,1,1))
    n=size(A)[1]
    for ii in 2:n
        for jj in 1:(ii-1) ## Just over one triangle
            if compV==false
                A=Rotate(A,ii,jj)
            else
                A,V=Rotate(A,ii,jj;computeV=true,V=V);
            end
        end
    end
    if compV==false
        return A
    else
        return A,V
    end
end

function jacobisweep(A,iters)
    n=size(A)[1]
    V = Matrix(Float64.(I(n)))
    Vout = Matrix(Float64.(I(n)))
    for i in 1:iters
        A,V = Sweep(A;compV=true,V)
        Vout = transpose(V)*Vout
    end
    return Matrix(A),Vout
end

#end
