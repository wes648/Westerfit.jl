"""
Okay this is the Hyperfine calculator. I'm going to build it for up to 9 spins,
    3 sets of 3 coupled

Procedure:
Based on F,S,I, initialize zeros matrix
    determine allowed Js from Δ(F,J,I)
        determine allowed Ns from Δ(J,S,N)
    make arrays of Jlvls, and Nlvsls
        construct Jlvls such that it's the same length as Nlvls
    sum(2 .* Nlvls .+ 1) for Fmat dims
for n in Nlvls
    for nb in n:Nlvls
        Fmat[dims:dims] = <FIJSNK|Hq|FIJ'S'N'K'> + <JSNK|Ho|J'S'N'K'>
    end
end
Fmat = Hermitian(Fmat)

adding spins:
For given F
for i in 1:length(Is)
    determine allowed Fn
    use I = sum(Is[end-i-1,:])
    allowed Fn from Δ(F,I,Fn)

"""

module hypcalc

using LinearAlgebra
using HalfIntegers
include("./WIGXJPF.jl")
using .WIGXJPF

function Δlist(J,I)
    max = J+I
    min = abs(J-I)
    return collect(min:max)
end

function hqelem(χtens,F,I,J,S,N,K,q,Jb,Nb,Kb)
    out = 0.0
    out += 1/(wig3j(I,2,I,-I,0,I))*sqrt((2*Jb+1)*(2*J+1))
    out *= wig6j(I,Jb,F,J,I,2)*sqrt((2*Nb+1)*(2*N+1))
    out *= ϵtns[k+1,q+3]*wig3j(Nb,2,N,-Kb,q,K)
    out *= (-1)^(2*J+I+F+2*Nb+S-Kb)
    return out
end

function hqmat(χtens,F,I,J,S,SRmats)
    js = Δlist(F,I)
    nlist = []
    jlist = []
    for i in length(js)
        list = Δlist(js[i],S)
        nlist = vcat(nlist,list)
        longj = fill(js[i],length(list))
        jlist = vcat(jlist,longj)
    end
    fdim = sum(2.0 .* nlist .+ 1.0)
    fmat = zeros(Float64,fdim,fdim)

end
#function HypCalc(Fmax,Is,S,jmats,paramtensors)
#    if typeof(sum(Is)+S)=Int64
#        Fmin = 0
#    else
#        Fmin = 1//2
#    end
#end

end
