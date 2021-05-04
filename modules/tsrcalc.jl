"""
This is TSRCalc!
This contains the actual novelty of the westersim project. It has the
spin-torsion-rotation hamiltonain and operator.
"""

module tsr
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearAlgebra.LAPACK
using Profile
using Distributed
using DelimitedFiles
using Printf
include("./rotcalc.jl")
using .rotcalc
include("./spicalc.jl")
using .spicalc
include("./qnassign.jl")
using .qn

S = 0.5
vtmax = 0
vtcalc = 8
vthardcap = 10
Nmax = 15



function vectrunator(vector::Array{Float64},N::Float64,J::Float64)::Array{Float64}
    if N>J
        shortvec = vector[(vtcalc+1)*2*N:end]
        shortvec = vcat(0.0,vector,0.0)
    else
        shortvec = vector[1:(vtcalc+1)*(2*N+2)]
    end
    shortvec = normalize!(shortvec)
    return shortvec
end
function Indexer(J,N,τ)
    oldjs = collect(0.5:(J-1))
    index = sum(4 .* oldjs .+ 2)
    index += (N-J+1/2)(2N-1) + N+τ+1
    return index
end
function findnearest(a::Array{Float64},x::Float64)
    length(a) > 0 || return 0:-1
    r = searchsorted(a,x)
    length(r) > 0 && return r
    last(r) < 1 && return searchsorted(a,a[first(r)])
    first(r) > length(a) && return searchsorted(a,a[last(r)])
    x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
    x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
    return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end
function closest_index(a::Array{Float64,1}, x::Float64)::Int64
    ibest = 1
    dxbest = abs(a[ibest]-x)
    for i in eachindex(a)
        dx = abs(a[i]-x)
        if dx < dxbest
            dxbest = dx
            ibest = i
        end
    end
    return ibest
end

########################################################
### Here is the coupled-torsional protion of H_{tot} ###
########################################################
function Ht0K0m0N(tvals::Array{Float64,2},tvecs::Array{Float64,3},tsrprm,
    thet::Float64,J::Float64,K::Float64,vt::Int64,sigma::Float64)::Float64
    out = zeros(Float64,2*vthardcap+1)
    #ms = collect(Float64,0:vtcalc)
    #@. out = eta*(3.0*ms+sigma)*K*thet*tvecs[ms,vt+1,Int64(K)+Nmax+1]^2
    #println("J=$J,K=$K,p=$(K+Nmax+1)")
    for m in 0:vtcalc#2*vthardcap
        mp = m+1
        out[mp] = tsrprm[1]*(3.0*m+sigma)*K*thet
        #the vt+1 here is to adjust for Julia's array indexing
        out[mp] *= tvecs[mp,vt+1,Int64(K)+Nmax+1]^2
        #out[mp] += tvals[mp,Int64(K)+Nmax+1]
    end
    out = sum(out) + tvals[vt+1,Int64(K)+Nmax+1]
    return out
end
function Ht0K0m1N(jphi::Float64,tsrprm,J::Float64,N::Float64,K::Float64,vt::Int64,
    sigma::Float64,tvecs)::Float64
    out = zeros(Float64,2*vthardcap+1)
    for m in -vthardcap:vthardcap
        mp = m+vthardcap+1
        out[mp] = tsrprm[1]*(3.0*m+sigma)
        out[mp] *= sqrt(N^2 - K^2)*jphi
        out[mp] *= tvecs[mp,vt+1,Int64(K)+Nmax+1]*tvecs[mp,vt+1,Int64(K)+Nmax+1]
    end
    out = sum(out)
    return out
end
#The below functions are only for higher order terms that I can turn on later
#function Ht1K0m0N(J,N,vt,sigma)
#function Ht2K0m0N(J,N,vt,sigma)
#function Ht0K1m0N(J,N,vt,sigma)
function Ht1Kov0N(J::Float64,N::Float64,K::Float64,vt::Int64,vb::Int64,
    sigma::Float64,tvecs::Array{Float64,3},rprms)::Float64
    Nt2 = N*(N+1)
    out = zeros(Float64,2*vthardcap+1)
    for m in -vthardcap:vthardcap
        mp = m+vthardcap+1
        out[mp]  = rprms[4]*sqrt(Nt2-K*(K-1.0))*(K-0.5)
        out[mp] *= tvecs[mp,vt+1,Int64(K)+Nmax+1]*tvecs[mp,vb+1,Int64(K)+Nmax]
    end
    out = sum(out)
    return out
end
function Ht2Kov0N(J::Float64,N::Float64,K::Float64,vt::Int64,vb::Int64,sigma::Float64,
    tvecs::Array{Float64,3},rprms)::Float64
    Nt2 = N*(N+1)
    out = zeros(Float64,2*vthardcap+1)
    for m in -vthardcap:vthardcap
        mp = m+vthardcap+1
        out[mp]  = 0.25*(rprms[2]-rprms[3]) + rprms[8]*Nt2 + 0.5*rprms[9]*(K^2 + (K-2.0)^2)
        out[mp] *= tvecs[mp,vb+1,Int64(K)+Nmax-1]
        out[mp] *= sqrt((Nt2-K*(K-1.0))*(Nt2-(K-1.0)*(K-2.0)))*tvecs[mp,vt+1,Int64(K)+Nmax+1]
    end
    out = sum(out)
    return out
end
#function Ht0K2m0N(J,N,vt,sigma)
#function Ht1K2m0N(J,N,vt,sigma)
#function Ht2K2m0N(J,N,vt,sigma)
#function Ht1K0m1N(J,N,vt,sigma)
#function Ht2K0m1N(J,N,vt,sigma)
#function Ht0K1m1N(J,N,vt,sigma)
#function Ht1K1m1N(J,N,vt,sigma)
#function Ht2K1m1N(J,N,vt,sigma)
#function Ht0K2m1N(J,N,vt,sigma)
#function Ht1K2m1N(J,N,vt,sigma)
#function Ht2K2m1N(J,N,vt,sigma)

function Hspitor0N(J::Float64,N::Float64,NI::Int64,v::Int64,sigma::Float64,
    thet::Float64,tvals::Array{Float64,2},tvecs::Array{Float64,3})::Array{Float64,2}
    #thet = theta(J,N)
    #karray = collect(Int64,-N:N)
    elems = zeros(Float64,2*NI+1)
    for k in -N:N
        elems[Int64(k)+NI+1] = Ht0K0m0N(tvals,tvecs,thet,J,N,k,v,sigma)
    end
    out = diagm(0=>elems)
    return out
end
function Hspitor0Nv(J::Float64,N::Float64,NI::Int64,vk::Int64,vb::Int64,sigma::Float64,
    tvecs::Array{Float64,3},rprms)::Array{Float64,2}
    if N==0
        out = [0.0]
    else
        o1elems = zeros(Float64,2*NI)
        o2elems = zeros(Float64,2*NI-1)
        o1elems[1] = Ht1Kov0N(J,N,-N+1,vk,vb,sigma,tvecs,rprms)
        for k in (-N+2):N
            o1elems[Int64(k)+NI] = Ht1Kov0N(J,N,k,vk,vb,sigma,tvecs,rprms)
            o2elems[Int64(k)+NI-1] = Ht2Kov0N(J,N,k,vk,vb,sigma,tvecs,rprms)
        end
        out = diagm(1=>o1elems,2=>o2elems)
        out = Symmetric(out)
    end
    return out
end
function Hspitor1N(J::Float64,N::Float64,NI::Int64,v,sigma::Float64,tvecs::Array{Float64,3}
    ,tsrprm)::Array{Float64,2}#this one is also fixed
    jphi = -1.0/(J+0.5)
    #Nf = convert(Float64,N)
    karray = collect(Float64,-N:N)
    elems = zeros(Float64,2*NI-1)
    for k in (-N+1):(N-1)
        elems[Int64(k)+NI] = Ht0K0m1N(jphi,tsrprm,J,N,k,v,sigma,tvecs)
    end
    tsmat = hcat(zeros(length(karray)-2),diagm(0=>elems),zeros(length(karray)-2))
    return tsmat
end

function Htot(J::Float64,Nf::Float64,Ni::Int64,sigma::Float64,
    tvals::Array{Float64,2},tvecs::Array{Float64,3},rprms::Array{Float64,1},
    srprms::Array{Float64,1},tsrprms)::Array{Float64,2}
#function Htot(J::Float64,Nf::Float64,Ni::Int64,rprms::Array{Float64,2},
#sigma::Float64,tvals::Array{Float64,2},tvecs::Array{Float64,3})::Symmetric{Float64,Array{Float64,2}}
    Nud = 2*Ni+1
    #println("N=$Ni")
    #println(rotcalc.Hrot(Nf,rprms))
    Nld = 2*Ni-1
    thet = spicalc.theta(J,Nf,S)
    Hsrumat = rotcalc.Hrot(Nf,rprms) + spicalc.Hspi0N(J,Nf,S,srprms)
    Hsrlmat = rotcalc.Hrot(Nf-1.0,rprms) + spicalc.Hspi0N(J,Nf-1.0,S,srprms)
    Hsrlmat = Hsrlmat[1:Nld,1:Nld]
    URsrmat = spicalc.Hspi1N(J,Nf,S,srprms)
    Numat = Array{Float64}(undef,0,(Nud)*(vtcalc+1))
    Nlmat = Array{Float64}(undef,0,(Nld)*(vtcalc+1))
    URmat = Array{Float64}(undef,0,(Nud)*(vtcalc+1))
    #This block here appears to be the rate limiting step
    for vt in 0:vtcalc
        vturow = Array{Float64}(undef,Nud,0)
        vtlrow = Array{Float64}(undef,Nld,0)
        urvtrow = Array{Float64}(undef,Nld,0)
        for vi in 0:vtcalc
            if vi==vt
                nutemp = Hsrumat + Hspitor0N(J,Nf,Ni,vt,sigma,thet,tvals,tvecs)
                nltemp = Hsrlmat + Hspitor0N(J,Nf-1.0,Ni-1,vt,sigma,thet,tvals,tvecs)[1:Nld,1:Nld]
                urtemp = URsrmat + Hspitor1N(J,Nf,Ni,vt,sigma,tvecs,tsrprms)
            elseif vt>vi
                nutemp = Hspitor0Nv(J,Nf,Ni,vt,vi,sigma,tvecs,rprms)
                nltemp = Hspitor0Nv(J,Nf-1.0,Ni-1,vt,vi,sigma,tvecs,rprms)[1:Nld,1:Nld]
                urtemp = zeros(Nld,Nud)
            else
                nutemp = zeros(Float64,Nud,Nud)
                nltemp = zeros(Float64,Nld,Nld)
                urtemp = zeros(Nld,Nud)
            end
            vturow = hcat(vturow,nutemp)
            vtlrow = hcat(vtlrow,nltemp)
            urvtrow = hcat(urvtrow,urtemp)
        end
        Numat = vcat(Numat,vturow)
        Nlmat = vcat(Nlmat,vtlrow)
        URmat = vcat(URmat,urvtrow)
    end
    LLmat = transpose(URmat)
    Htotmat = [Nlmat URmat; LLmat Numat]
    #Htotmat = Symmetric(Htotmat)
    return Htotmat
end
function Htotf0(J::Float64,sigma::Float64,tvals::Array{Float64,2},
    tvecs::Array{Float64,3},rprms::Array{Float64,1},
    srprms::Array{Float64,1},tsrprms)::Array{Float64,2}
    thet = 0.5#theta(0.5,1.0)
    Hsrumat = rotcalc.Hrot(1.0,rprms) + spicalc.Hspi0N(J,1.0,S,srprms)
    Hsrlmat = zeros(Float64,1)
    URsrmat = spicalc.Hspi1N(J,1.0,S,srprms)   #wildly enough this line specifically is super slow
    Numat = Array{Float64}(undef,0,3*(vtcalc+1))
    Nlmat = Array{Float64}(undef,0,1*(vtcalc+1))
    URmat = Array{Float64}(undef,0,3*(vtcalc+1))
    for vt in 0:vtcalc
        vturow = Array{Float64}(undef,3,0)
        vtlrow = Array{Float64}(undef,1,0)
        urvtrow = Array{Float64}(undef,1,0)
        for vi in 0:vtcalc
            if vi==vt
                nutemp = Hsrumat + Hspitor0N(J,1.0,1,vt,sigma,thet,tvals,tvecs)
                nltemp = Hsrlmat + Hspitor0N(J,0.0,0,vt,sigma,thet,tvals,tvecs)
                urtemp = URsrmat + Hspitor1N(J,1.0,1,vt,sigma,tvecs,tsrprms)
            else
                nutemp = Hspitor0Nv(J,1.0,1,vt,vi,sigma,tvecs,rprms)
                nltemp = zeros(Float64,1) #Hspitor0Nv(J,0,vt,vi,sigma,tvecs)
                urtemp = zeros(1,3)
            end
            vturow = hcat(vturow,nutemp)
            vtlrow = hcat(vtlrow,nltemp)
            urvtrow = hcat(urvtrow,urtemp)
        end
        Numat = vcat(Numat,vturow)
        Nlmat = vcat(Nlmat,vtlrow)
        URmat = vcat(URmat,urvtrow)
    end
    LLmat = transpose(URmat)
    Htotmat = [Nlmat URmat; LLmat Numat]
    #Htotmat = Symmetric(Htotmat)
    return Htotmat
end

function TSRDiag(Nf::Float64,Ni::Int64,sigma::Float64,
    tvals::Array{Float64,2},tvecs::Array{Float64,3},rprms,srprms,tsrprms)#::Array{Float64,2}
    J = Nf-0.5
    #println("J=$J, matrix construction")
    Hmat = Htot(J,Nf,Ni,sigma,tvals,tvecs,rprms,srprms,tsrprms)
    #println("eigens then qn assigns")
    #@time Hmat = eigen!(Hmat)
    Hvals, Hvecs = LAPACK.syev!('V','U', Hmat)
    #output = qnlocalglobal(N, rprms,sprms,sigma, Hmat.values, Hmat.vectors, tvals, tvecs)
    #@time output = qn.suboverlap(Nf,vtcalc,vtmax,Hmat.values,Hmat.vectors)
    Hvals, Hvecs, QNs = qn.suboverlap(Nmax,Nf,vtcalc,vtmax,Hvals,Hvecs,sigma)
    return Hvals, Hvecs, QNs
end

function TSRDiagf0(sigma::Float64,tvals::Array{Float64,2},tvecs::Array{Float64,3},
    rprms,srprms,tsrprms)#::Array{Float64,2}
    #J = 0.5
    Hmat = Htotf0(0.5,sigma,tvals,tvecs,rprms,srprms,tsrprms)
    #Hmat = eigen!(Hmat)
    Hvals, Hvecs = LAPACK.syev!('V','U', Hmat)
    #output = qnlocalglobal(1, rprms,sprms,sigma, Hmat.values, Hmat.vectors, tvals, tvecs)
    #output = qn.suboverlap(1,vtcalc,vtmax,Hmat.values,Hmat.vectors)
    Hvals, Hvecs, QNs = qn.suboverlap(Nmax,1,vtcalc,vtmax,Hvals,Hvecs,sigma)
    return Hvals, Hvecs, QNs
end

function RotSpiTorCalc(Nmax::Int64,sigma::Float64,tvals::Array{Float64,2},tvecs::Array{Float64,3},
    rprms,srprms,tsrprms)#::Array{Float64,2}
    marray = [-vthardcap:vthardcap]
    niarray = collect(Int64,2:Nmax)
    nfarray = collect(Float64,2:Nmax)
    egys, wvfn, qunu = TSRDiagf0(sigma,tvals,tvecs,rprms,srprms,tsrprms)
    for i in 1:Nmax-1
        #println("Nl=$i")
        newe, neww, newq = TSRDiag(nfarray[i],niarray[i],sigma,tvals,tvecs,rprms,srprms,tsrprms)
        egys = vcat(egys,newe)
        wvfn = hcat(wvfn,neww)
        qunu = vcat(qunu,newq)
    end
    qunu = Int64.(qunu)
    return egys, wvfn, qunu
end



end
