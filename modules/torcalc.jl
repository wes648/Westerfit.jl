"""
This sub-program calculates the torsional energy levels and vectors before the
    rotational calculations. It is called TorCalc
Basis set: Ψ=|Kmσ⟩

Potential improvements:
    n-fold rotor treatement
    flexible parameter definitions
"""

module torcalc
using LinearAlgebra
vthardcap = 10

function vecvalpair(m::Array{Int64})::Array{Int64,1}
#This function provides the definitional pairing of the torsional eigenvalues to the eigenvector elements
	out = @. vthardcap+1 - div((m-(m%2))*((-1)^(m%2)),2)
#    out = @. 11 + out
    return out
end

#Matrix Elements by Δm
#torparams = [F; rho; V3; V6; fv3; fv6]
function Htorm0(K::Float64,m::Array{Float64,1},σ::Float64,tprms)::Array{Float64,1}
#	out = @.        F*(3.0*m+σ-rho*K)^2 + V3*0.5 + V6*0.5 + 0.0*fv3*(3.0*m+σ)^2 + 0.0*fv6*(3.0*m+σ)^2
	out = @. tprms[1]*(3.0*m+σ-tprms[2]*K)^2 + tprms[3]*0.5 + tprms[4]*0.5
	out = @. out + 0.0*tprms[5]*(3.0*m+σ)^2 + 0.0*tprms[6]*(3.0*m+σ)^2
end
function Htorm1(K::Float64,m::Array{Float64,1},σ::Float64,tprms)::Array{Float64,1}
    out = @. -tprms[3]*0.25 + 0.0*tprms[5]*(3.0*m+σ)^2
end
function Htorm2(K::Float64,m::Array{Float64,1},σ::Float64,tprms)::Array{Float64,1}
    out = @. -tprms[4]*0.25 + 0.0*tprms[6]*(3.0*m+σ)^2
end


function Htors(K::Float64,σ::Float64,tprms)#::Tuple{Array{Float64,2},Array{Float64,3}}
#this function builds the torsional matrix for a given K & σ
    marray = collect(Float64,-vthardcap:vthardcap)
    ondiags = Htorm0(K,marray,σ,tprms)
    of1diag = Htorm1(K,marray[1:end-1],σ,tprms)
    of2diag = Htorm2(K,marray[1:end-2],σ,tprms)
    Tormat = LinearAlgebra.diagm(0=>ondiags, 1=>of1diag, 2=>of2diag)
    #Tormat = Hermitian(Tormat)
    #Tormat = eigen!(Tormat)
	vals, vecs = LAPACK.syev!('V','U', Tormat)
	#if K==0
	#	(Tormat.values/29979.2458)
	#else
	#	nothing
	#end
#    return Tormat.values,Tormat.vectors
    return vals, vecs
end
function tormatch(vects::Array{Float64,2})::Array{Float64,2}
#This reorders the
    test = collect(Int64,1:2*vthardcap+1)
    newindices = vecvalpair(test)
    vects[test,:] = vects[newindices,:]
    return vects
end
#function TorCalc(parameters::Array{Float64,2}, settings::Array{Float64},
#Nm::Float64,σ::Float64)::Tuple{Array{Float64,2},Array{Float64,3}}
function TorCalc(parameters::Array{Float64,1}, settings, Nm::Int64,
	σ::Float64)::Tuple{Array{Float64,2},Array{Float64,3}}
    vthardcap = settings[1]
    #Okay so I need to reorder the eigenvectors as well.
    #I will out put a 3D array of eigenvectors and a 2D of values.
	#dimenions are m,vt,k
    #This will not be truncated. Truncation will occur in what is called in the next stage
    #marray = collect(1:21)
	Nf = Float64(Nm)
	torvals = zeros(Float64,2*vthardcap+1)
	torvecs = zeros(Float64,2*vthardcap+1,2*vthardcap+1)
	torvals,torvecs = Htors(-Nf,σ,parameters)
    torvecs = tormatch(torvecs)#tores,torvs)#,marray)
	karray = collect(Float64,-Nf+1.0:Nf)
    for k in karray
        newvals,newvecs = Htors(k,σ,parameters)
        newvecs = tormatch(newvecs)#tores,torvs)#,marray)
        torvals = hcat(torvals,newvals)
        torvecs = cat(torvecs,newvecs,dims=3)
    end
    return torvals,torvecs
end

end
