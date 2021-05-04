"""
This will serve as a test form for TorSpiSymCalc with a hardcoded Hamiltonian
A generalized from of SymCalc currently exists. This should make the creation
of WesterFit easier. Or at least I hope it will

There are still some errors in the energy levels as compared to SPFIT & BELGI
    Not sure what the source of error is. Should break down by individual parameters

I should remove as many type changes as I can to speed up the code

Current QN approach is pretty fun
    I'm taking the overlap integral for each vt state and finding the largest magnitue element
    This is being used to assign which N and vt the state is.
    I'm still just assigning τ by sorting. Should consider similarity factors

For now I need to work on the simulation stage. It could be cleaned up massively
    Unfortunately, RAM36's check everything approach is what I'm going to need to adopt here
    At least I'll finally get myself to learn the details of line intensities

Spin-torsion-rotation is pretty great
"""

#println("")
#using GenericLinearAlgebra
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearAlgebra.LAPACK
using Profile
using Distributed
using DelimitedFiles
using Printf
include("./modules/torcalc.jl")
using .torcalc
include("./modules/rotcalc.jl")
using .rotcalc
include("./modules/spicalc.jl")
using .spicalc
include("./modules/tsrcalc.jl")
using .tsr
include("./modules/qnassign.jl")
using .qn
include("./modules/tracalc.jl")
using .tracalc

#c = 29979.2458 #MHz to cm-1 converter

#Etors = zeros(5,3)
## Below are some general purpose functions ##
#function f(x::Array{Float64},y::Array{Float64})::Float64

#function vecvalpair(m::Array{Int64})
#This function provides the definitional pairing of the torsional eigenvalues to the eigenvector elements
#    out = @. vthardcap+1 - div((m-(m%2))*((-1)^(m%2)),2)
#    out = @. 11 + out
#    return out
#end
function fulleigencalc(rprms::Array{Float64,1},sprms::Array{Float64,1},torparams::Array{Float64,1},
    tsrparams,Nmax::Int64,torsets)#::Array{Float64,3}
    if useTor==true
        σ0val,σ0vec = torcalc.TorCalc(torparams, torsets, Nmax, 0.0)
        σ1val,σ1vec = torcalc.TorCalc(torparams, torsets, Nmax, 1.0)
        s0val,s0vec,s0qns = tsr.RotSpiTorCalc(Nmax,0.0,σ0val,σ0vec,rprms,sprms,tsrparams)
        s1val,s1vec,s1qns = tsr.RotSpiTorCalc(Nmax,1.0,σ1val,σ1vec,rprms,sprms,tsrparams)
        EgyWriter(s0val,s0qns,0)
        EgyWriter(s1val,s1qns,1)
    else
        σ0val = zeros(2*vthardcap+1,2*Nmax+1)
        σ0vec = zeros(2*vthardcap+1,2*vthardcap+1,2*Nmax+1)
        σ1val = zeros(size(σ0val))
        σ1vec = zeros(size(σ0vec))
        s0val,s0vec,s0qns=RotSpiTorCalc(Nmax,sprms,0.0,rprms,srprms,tsrparams)
        sig1ful=zeros(size(sig0ful))
        EgyWriter(s0val,s0qns,0)
    end
    finaleigs = cat(s0val,s1val, dims=3)
    finalvecs = cat(s0vec,s1vec, dims=4)
    finalquns = cat(s0qns,s1qns, dims=4)
    torvecs = cat(σ0vec,σ1vec, dims=4)
    return finaleigs, finalvecs, finalquns, torvecs
end

function EgyWriter(energies, qunus,sigma)
    c = 29979.2458
    out = fill("0",size(energies)[1])
    for i in 1:size(energies)[1]
        energy = energies[i]/c
        part = lpad(@sprintf("%0.4f", energy), 15)
        part = string(part,lpad(@sprintf("%0.1f", qunus[i,1,1]/2), 7))
        out[i] = string(part, lpad(qunus[i,2,1],4),lpad(qunus[i,3,1],4), lpad(qunus[i,4,1],4),
        lpad(qunus[i,5,1],4), lpad(qunus[i,6,1],4))
#        lpad(qunus[i,5,1],4), lpad(sigma,4))
    end
    if sigma==0
        io = open("Astates_$molnam.egy", "w") do io
            for i in out
                println(io, i)
            end
        end
    elseif sigma==1
        io = open("Estates_$molnam.egy", "w") do io
            for i in out
                println(io, i)
            end
        end
    else
        println("Sorry, not ready for this σ value yet")
    end
end


#function westersim()
useSpi=true
useTor=true
useNatUnits=false
S = 0.5
A = 0.0
B = 0.0
C = 0.0
Dab = 0.0
DelNK = 0.0
DelK = 0.0
DelN = 0.0
delN = 0.0
delK = 0.0
F = 0.0
V3 = 0.0
V6 = 0.0
rho = 0.0
fv3 = 0.0
fv6 = 0.0
ezz = 0.0
exx = 0.0
eyy = 0.0
exz = 0.0
ezx = 0.0
DelsK = 0.0
DelsN = 0.0
DelsNK = 0.0
DelsKN = 0.0
delsN = 0.0
delsK = 0.0
eta = 0.0
μa = 1.0
μb = 1.0
vtmax = 0
vtcalc = 8
vthardcap = 10
Nmax = 15
#allowed = ["S", "A", "B", "C", "Dab", "DelNK", "DelK", "DelN", "delN",
#"delK", "F", "V3", "V6", "rho", "fv3", "fv6", "ezz", "exx", "eyy",
#"exz", "ezx", "DelsK", "DelsN", "DelsNK", "DelsKN", "delsN", "delsK",
#"eta", "μa", "μb", "vtmax", "vtcalc", "Nmax", "vthardcap"]
#You may have noticed the absurd number of νt limits.
#Here is what they mean in the order they are used in the code:
#νthardcap is the value used in the 1st diagonalization stage.
#   BELGI uses 10. It has a silly name because it's probably not worth
#   futzing with unless you are calculating very high torsional states
#vtcalc is the highest vt in the 2nd diagonalization stage.
#   The truncation is the improve calculation speeds.
#   Should probably be at least 3 higher than what you want in the output
#vtmax is the highest vt that is output by the program.
println("Westersim!")
print("Enter molecule name: ")
println("")
#molnam = readline()
molnam = "spitorrot"
include("$molnam.inp")
println(molnam)
println("Nmax=$Nmax")
ρ=rho
Dab = 0.5*Dab
#    global vthardcap
#    global vtcalc
#    global vtmax
#    global Nmax
#    global S
#    global molnam
    #if useNatUnits==true
    #    F = 29979.2458*F
    #    V3 = 29979.2458*V3
    #else
    #    nothing
    #end
    #if useSpi==true
    #    ao = -(eaa+ecc+ebb)/3.0
    #    a = -(2.0*eaa-ecc-ebb)/6.0
    #    d = -(ezx + exz)*0.5
    #    b = -(ebb-ecc)*0.5
    #else
    #    ao = 0.0
    #    a = 0.0
    #    d = 0.0
    #    b = 0.0
    #end
function westersim()
    #if useTor==false
    #    vtmax = 0
    #    vtcalc = 0
    #    ρ = 0.0
    #    F = 0.0
    #    V3 = 0.0
    #else
    #    thetaRAM = 0.5*atan(2*Dab/(A-B))
    #    μa = μa*cos(thetaRAM) + μb*sin(thetaRAM)
    #    μb = -μa*sin(thetaRAM) + μb*cos(thetaRAM)
    #end
    torparams = [F; rho; V3; V6; fv3; fv6]
    torsets = [vthardcap, vtcalc, vtmax]
    rotparams = [A; B; C; Dab; DelN; DelNK; DelK; delN; delK]
    ao = -(ezz+eyy+exx)/3.0
    a = -(2.0*ezz-eyy-exx)/6.0
    d = -(ezx + exz)*0.5
    b = -(exx-eyy)*0.5
    spinrotparams = [ao; a; b; 0.0; d; 0.0; DelsN;DelsK; DelsNK; DelsKN; delsN; delsK]
    strp = [eta]
    #println(rotparams)
    #println(torparams)
    #println(spinrotparams)
    println("Calculating Energy Levels")
    @time fegy, fvec, fqun, ftvcs = fulleigencalc(rotparams,spinrotparams,torparams,strp,Nmax, torsets)
    #println("Calculating A State Transitions, then E states")
    #println(size(fvec))
    println("Calculating A State Transitions")
    @time Afreqs,Aqunus = tracalc.TraCalcCs(Nmax,S,vtcalc,fegy[:,:,1],fvec[:,:,1,1],
    ftvcs[:,:,:,1],fqun[:,:,:,1],0,0)
    if useTor==true
        println("Calculating E State Transitions")
        @time Efreqs,Equnus = tracalc.TraCalcCs(Nmax,S,vtcalc,fegy[:,:,2],fvec[:,:,1,2],
        ftvcs[:,:,:,2],fqun[:,:,:,2],1,0)
        ffreqs = vcat(Afreqs,Efreqs)
        fquns = vcat(Aqunus, Equnus)
        println("Writing Transitions to file")
        @time tracalc.TraWriter(molnam,ffreqs, fquns)
    else
        println("Writing Transitions to file")
        @time tracalc.TraWriter(molnam, Afreqs, Aqunus)
    end
end

#There is still plenty of optimizations to do here
#also need to figure out what's up with σ=±1

useTor=true
#println(Nmax)

@time westersim()
println()
println("A miracle has come to pass. The end of the code was reached!")
#println()
