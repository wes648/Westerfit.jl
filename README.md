# westerfit
A new program for the simulating and fitting of molecular rotational spectra for $C_s$ molecules with one internal rotor and one spin source.

The westerfit dev team is back! Questions can be directed to wes.harper2@gmail.com

**WARNING** This code has not yet been published! You are more than welcome to download it & fiddle around with it but if you need to publish the results, please reach out to the authors first! Once the paper is out, we will remove this warning and you'll be able to publish results without contacting us first.

## Quickstart

The following is intended to be a quick reference on the input file structure & program usage. A more complete manual will be constructed once the paper is published.
The westerfit input file is divided into three sections each with a designated header. The very first line is a title card and the sections are the Control Settings (%CNTRLS), Second Order Parameters (%2NDORDER), and the additional Operators & Parameters (%PARAMS).

### Control Settings
Below are control settings, their meanings, and default values::Type

**apology** (true::Bool): Prints a sincere if awkward apology for the name of the program

**RUNmode** (ESF::String): Dictates how the program is run by finding characters in string. E calculates & prints energy levels, S calculates & prints transitions, and F performs the optimization. ES will run through the same calculations as S but without printing energy levels. The Fit will run first and then put the new parameters into the ES calculation.

**overwrite** (true::Bool): Whether or not to overwrite the input file with the new values from the fit. Regardless, the input file will be copied to a back up.

**S** (0::Float64): Spin value of the internal spin source. This can be either electron or nuclear spin

**NFOLD** (0::Int): Symmetry-fold of the internal rotor. 0 disables most of the torsional code (I think)

**mcalc** (8::Int): Determines the size of the torsional basis as 2mcalc+1 for A & E symmetry states and 2mcalc+2 for B symmetry. BELGI uses 10 for the first stage and the equivalent of 4 for the second stage. I've found 5 to work reasonably well without being repulsively slow.

**vtmax** (0::Int): Largest torsional state output by the code. I *think* this also impacts what torsional states are used in the fitter. Best to keep with just ground state for now. 

**Jmax** (0.0::Float64): Largest J value in the simulation. The code will automatically adjust if S is a half integer and print a warning.

**νmin** (0.0::Float64): Lowest frequency in the simulation in GHz

**νmax** (40.::Float64): Highest frequency in the simulation in GHz. Make sure this is larger than νmin or the code will crash

**INTthres** (0.0::Float64): Minimum intensity value included in the file

**TK** (8.0::Float64): Temperature of the simulation, not sure this is actually used because a certain author finds thermodynamics abhorrent

**λlm0** (0.0001::Float64): Scale-factor used to determine the inital Levenberg-Marquardt Parameter. This gets multiplied by a function of the rms to determine the LBMQ Parameter used in a given step.

**turducken** (1::Int): Number of Levenberg-Marquardt steps calculated on single step before recalculated the Jacobian. Doesn't seem worth it given the current performance of the Jacobian but can give an occasional performance boost

**assign** (expect::String): Determines how the quantum numbers are assigned after diagonalization. The default, expect, uses the expectation values of $m$ & then $N$ followed by an energetic sorting to assign $K$. RAM36 uses the contributions of different blocks of the eigenvectors to provide a spin-expanded version of the assigner in RAM36. expectk is similar to expect but uses the expectation values of $K$ as well and doesn't seem to work. Lastly, eeo does the expection values of $m$ and $N$ followed by eigenvector analysis to assign $K$. This does the best job of reproducing SPFIT's DIAG=3. Generally expect is recommended but RAM36 works very nicely for single perturbations (spin or torsion). I'm personally fond of the theory in eeo but find its performance lacking.


### Second Order Parameters
The second order Hamiltonian of westerfit is hardcoded. You can comment out any lines with # but the number of lines must remain fixed in this section. All the parameters default to zero and there are some internal transformations that occur upon initializing the code.

$A$, $B$, and $C$ are the rotational constants corresponding to the $z$, $x$, and $y$ molecule fixed axes respectively. 
$D_{ab}$ is the off-diagonal rotational constant for the Rho Axis Method. 
When it is included, $A$ and $B$ refer to their RAM meanings not their PAM meanings. These are all in MHz

$F$ is the internal rotor constant. It is the reduced value as is used in RAM36. 
$\rho$ is the coupling term between the methyl rotation and the $z$-axis angular momentum. 
It engages in the Hamiltonian as $-2\rho FP_{\alpha}N_{z}$ and $(A-F\rho^{2})N_{z}^{2}$. 
The $Vn$ is the first term in the potential expansion as determined by NFOLD.

The $\epsilon$ terms are the spin-rotation coupling terms referring to axes in their subscripts. 
ϵxz is the average of the $\epsilon_{xz}$ and $\epsilon_{zx}$ terms as they do not have seperable matrix elements. 
These are transformed into spherical tensor notation upon start of the code.
These can also serve as nuclear spin-rotation coupling as they are mathematically equivalent.

The $\chi$ terms are the quadrupole terms, again referring to the axes in their subscript.
They are also transformed into spherical tensor form upon initializaiton.

Lastly is the spin-torsion coupling term $\eta$.

### Higher order operators
These are manual coded in operators that are implemented as the anti-commutator of what the user codes in.
These lines can also be commented out but do not remove the lines opening with %.
The first column is a name string for the operators.
Here are some unicode characters for easier name: Δ, δ, Φ, ϕ, ϵ, χ, ρ, η, μ.
The second & third columns are Float64 and are the parameter value and step scale factor, respectively. A scale factor of zero will keep that operator fixed during the optimization.
The next 8 columns are Int referring to the powers of the various operators as described in the line beginign that block of the input file.
There are no checks of symmetry like in RAM36 so go wild. 
One could also code in inverse powers but I'm not sure why one would. Let me know if you do and how it helped!
The last column is a stage. It is either 0 for intensites of 1 for Hamiltonian operators. Might expand that if I come up with better code structures.


## Installation
Currently, the installation of westerfit is a touch convoluted.
You will need to install [Julia](https://julialang.org/) and I recommend doing so through [juliaup](https://github.com/JuliaLang/juliaup). 
Add the SparseArrays, StaticArrays, and WIGXJPFjl packages to your Julia installation.
The WIGXJPFjl package is a wrapper for [WIGXJPF](http://fy.chalmers.se/subatom/wigxjpf/) and I don't know if the install script will work on Windows.
You may have to manually place the shared library in the packages deps directory and I ask that you let me know if that's the case.
You may also have better luck using WSL on Windows if you don't have access to a Linux machine.
The code also uses DelimitedFiles, LinearAlgebra, Printf, and Dates but I believe all of those are included in Base.
Lastly, I recommend making a bash script called `westerfit` somewhere in your path.
This script just needs the following two lines:
```
 #!/bin/bash
 julia -tX /path/to/westerfit/src/run.jl $1
```
Set X to be the number of threads you want westerfit to run on (more is better!) and fill in your full path to the code.
Now all you need to do to run westerfit is type `westerfit molnam` and it'll run on molnam.inp.
You can replace molnam with any string.
I recommend a helpful molecule name.
Enjoy!


## A remark of inspiration:

>I don't have any idea whether we will ever need spectroscopy again in 40 years. 
>As for spectroscopy, people who want to do spectroscopy seem to be born that way.... 
>As I describe it, it's just one step in the compulsive/obsessive spectrum before the stage where you're put in an institution for your entire life. 
>And they have that, and they just love this stuff, and you can't get them to stop talking. 
>They're only 25 years old, I mean, they should be drinking beer and doing whatever, they want to talk about this and they love it. 
>So, those people are going to continue to do this, in the backwater universities, not heavily funded, but they are going to preserve the skills, because they can't help themselves. 
>They're born that way.
>
>And they simply don't want to sell cars. 
>They don't want to broad-brush the future of American science. 
>They want to do something that has eight digit numbers and a theory and obs minus calcs, and you're born that way. 
>You can see it in them.
>
>-Jon T. Hougen

