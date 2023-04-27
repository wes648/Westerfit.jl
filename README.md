# westerfit
A new program for the simulating and fitting of molecular rotational spectra for $C_s$ molecules with one internal rotor and one spin source.

Description of current status:
Torsion-Rotation tested against BELGI & RAM36. See test/belgi and example/2ba
Hyper-fine tested against SPFTI. see test/spfit
SR undergoing testing. see test/spfit

The westerfit input file is divided into three sections each with a designated header. The very first line is a title card and the sections are the Control Settings (%CNTRLS), Second Order Parameters (%2NDORDER), and the additional Operators & Parameters (%PARAMS).

### Control Settings
Below are control settings, their meanings, and default values::Type


**apology**: Prints a sincere if awkward apology for the name of the program, true::Bool

**RUNmode**: Dictates how the program is run by finding characters in string. E calculates & prints energy levels, S calculates & prints transitions, and F performs the optimization. ES will run through the same calculations as S but without printing energy levels. Currently ES run before F but this will be reversed in a future update. ESF::String


**S**: Spin value of the internal spin source. This can be either electron or nuclear spin, 0::Float64

**NFOLD**: Symmetry-fold of the internal rotor. 0 disables most of the torsional code (I think),  0::Int

**mcalc**: Determines the size of the torsional basis as 2mcalc+1 for A & E symmetry states and 2mcalc+2 for B symmetry. BELGI uses 10 for the first stage and the equivalent of 4 for the second stage. 8::Int

**vtmax**: Largest torsional state output by the code. I *think* this also impacts what torsional states are used in the fitter. Best to keep with just ground state for now. 0::Int


**Jmax**: Largest J value in the simulation, 0.0::Float64

**νmin**: Lowest frequency in the simulation in GHz, 0.0::Float64

**νmax**: Highest frequency in the simulation in GHz, 40.::Float64

**INTthres**: Minimum intensity value included in the file, 0.0::Float64

**TK**: Temperature of the simulation, not sure this is actually used because a certain author finds thermodynamics abhorrent, 8.0::Float64


**λlm0**: Scale-factor used to determine the inital Levenberg-Marquardt Parameter. This gets multiplied by rms/(1+rms), 0.0001::Float64

**turducken**: Number of Levenberg-Marquardt steps calculated on single step before recalculated the Jacobian. Doesn't seem worth it given the current performance of the Jacobian, 1::Int


### Second Order Parameters
The second order Hamiltonian of westerfit is hardcoded. You can comment out any lines with # but the number of lines must remain fixed in this section. All the parameters default to zero and there are some internal transformations that occure upon initializing the code.

$A$, $B$, and $C$ are the rotational constants corresponding to the $z$, $x$, and $y$ molecule fixed axes respectively. 
$D_{ab}$ is the off-diagonal rotational constant for the Rho Axis Method. 
When it is included, $A$ and $B$ refer to their RAM meanings not their PAM meanings. These are all in MHz

$F$ is the internal rotor constant. It is the reduced value as is used in RAM36. 
$\rho$ is the coupling term between the methyl rotation and the $z$-axis angular momentum. 
It engages in the Hamiltonian as $-2\rho FP_{\alpha}N_{z}$ and $(A-F\rho^{2})N_{z}^{2}$. 
The $Vn$ is the first term in the potential expansion as determined by NFOLD.

The $\epsilon$ terms are the spin-rotation coupling terms referring to axes in their subscripts. 
$\epsilon$xz is the average of the xz and zx terms as they do not have seperable matrix elements. 
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
The second & third columns are Float64 and are the parameter value and step scale factor, respectively. Currently the step scale factor just acts as a binary of 0 meaning don't fit or non-zero meaning do fit. Will be expanded later.
The next 8 columns are Int referring to the powers of the various operators as described in the line beginign that block of the input file.
There are no checks of symmetry like in RAM36 so go wild. 
One could also code in inverse powers but I'm not sure why one would. Let me know if you do and how it helped!
The last column is a stage. It is either 0 for intensites of 1 for Hamiltonian operators. Might expand that if I come up with better code structures.



### A remark of inspiration:

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

-Jon T. Hougen

