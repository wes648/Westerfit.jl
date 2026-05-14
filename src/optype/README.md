Please note this readme is for the up coming update for the version of westerfit with reworked guts & inputs.

# westerfit
A new program for the simulating and fitting of molecular rotational spectra for $C_s$ molecules with one internal rotor and one spin source.
The westerfit package has been developed & is maintained by J.H. "Wes" Westerfield and Sophie E. Worthingon-Krisch.

The paper is available [here](https://doi.org/10.1016/j.jms.2024.111928) and the pre-print is available [here](https://dx.doi.org/10.2139/ssrn.4807560).

Please feel free to direct any questions about the program to westerfit@proton.me

# Quickstart

The following is intended to be a quick reference on the input file structure & program usage. 
A more complete manual will be constructed once the paper is published.
The new westerfit input file is based on the [TOML](https://toml.io/en/) file structure for both human & machine readability.
The westerfit input file is divided into 5 sections each with a designated header.
The order of these blocks does not matter.
The sections will be summarized here and described in detail below.
The first \[info\] is written by the program and contains the molecule name (molnam), time of last run, and other information.
The \[controls\] section handles run behavior for the program as well as defining certain global paramters like S and N_fold.
The \[second_order\] block contains the information for the hard coded Hamiltonian.
Since these parameters are very commonly called, they are written with especially optimized routines.
The \[user_def\] section defines additional Hamiltonian operators including distortion terms.
Lastly, the \[intensities\] section defines the dipole operators.

### Control Settings
Below are control settings, their meanings, and default values::Type

**apology** (true::Bool): Prints a sincere if awkward apology for the name of the program

**RUNmode** (ESF::String): Dictates how the program is run by finding characters in string. 
E calculates & prints energy levels, S calculates & prints transitions, and F performs the optimization. 
ES and S will run through the same calculations but S will only print transitions while ES will also prints energy levels. 
The Fit will run first and then put the new parameters into the ES calculation.

**stages** (0::Int): Determines the number of diagonalization stages.
1 stage handles everything in a single matrix per J σ pair.
It is only available for 0 or 1 tops and only recommended for 0 tops.
2 stages treats a pure torsional stage per σ and then handles each J σ pair in a combined torsion-rotation stage.
This is usable for 1 or more tops.
3 stages treats each individual top separately, then builds a top-top stage, and then finally builds a torsion-rotation stage.
This is recommended for 2 or more tops.
The default value of 0 will cause a crash so this must be defined by the user.

**overwrite** (true::Bool): Whether or not to overwrite the input file with the new values from the fit.
Regardless, the input file will be copied to a back up.

**NFOLD** ([0]::Vector{Int}): Symmetry-fold of the internal rotors.
For multiple rotors, they don't have to have the same symmetry folds.
For example, para-cresol (HO-Ph-CH3) would be [3;2].
[0] disables most of the torsional code (I think).

**S** (0::Float64): Spin value of the internal spin source. This can be either electron or nuclear spin

**Jmax** (0.0::Float64): Largest J value in the simulation.
The code will automatically adjust if S is a half integer and print a warning.

**mcalc** (10::Int): Determines the size of the torsional basis as 2mcalc+1 for the first stage in which the rotors appear.
For 1 stage mode, this is applied to the torsion-rotation stage.
For 2 stage mode, the top-top stage is of size (2mcalc+1)^(no. of tops).
For the 3 stage mode, the single top stage is of size 2mcalc+1. 
BELGI uses 10 for the first stage and the equivalent of 4 for the second stage. 

**vtcalc** (8::Int): Determines the number of torsional states perserved after diagonalization.
This is not used for 1 stage mode.
A total of vtcalc+1 torsional stages are used for the torsion-rotation stage for both 2 and 3 stages mode.
In the 3 stage mode, vtcalc+1 torsional stages are taken from each indivudal top to build a (vtcalc+1)^(no. of rotors) sized top-top matrix.

**vtmax** (0::Int): Largest torsional state output by the code.
I *think* this also impacts what torsional states are used in the fitter.
Best to keep with just ground state for now. vtmax=1 is necessary for 2-fold rotors. I NEED TO CHECK THE BEHAVIOR

**νmin** (0.0::Float64): Lowest frequency in the simulation in GHz. The first character is a *nu* not a v so be careful

**νmax** (40.::Float64): Highest frequency in the simulation in GHz. Make sure this is larger than νmin or the code will crash

**INTthres** (0.0::Float64): Minimum intensity value included in the file

**TK** (8.0::Float64): Temperature of the simulation in Kelvin

**λlm0** (0.0001::Float64): Scale-factor used to determine the inital Levenberg-Marquardt Parameter. 
This gets multiplied by a function of the rms to determine the LBMQ Parameter used in a given step.

**turducken** (1::Int): Number of Levenberg-Marquardt steps calculated on single step before recalculated the Jacobian. 
Doesn't seem worth it given the current performance of the Jacobian but can give an occasional performance boost

**goal** (1.0::Float64): The value of the weighted RMS that terminates the fit. Intended to prevent fitting past the experimental resolution.

**assign** (ram36::String): Determines how the quantum numbers are assigned after diagonalization.
Currently only **ram36** works but others may get (re-)implemented later.
<!-- The default, **expect**, uses the expectation values of $m$ & then $N$ followed by an energetic sorting to assign $K$. 
**RAM36** uses the contributions of different blocks of the eigenvectors to provide a spin-expanded version of the assigner in RAM36. 
**expectk** is similar to expect but uses the expectation values of $K$ as well and doesn't seem to work. 
Lastly, **eeo** does the expection values of $m$ and $N$ followed by eigenvector analysis to assign $K$. 
This does the best job of reproducing SPFIT's DIAG=3. Generally **expect** is recommended but **RAM36** works very nicely for single perturbations (spin or torsion). 
I'm personally fond of the theory in **eeo** but find its performance lacking. -->

### Second Order Parameters
The second order Hamiltonian of westerfit is hardcoded. 
All the parameters default to zero and there are some internal transformations that occur upon initializing the code.
Each line must be in the structure of: <br>
N = [val, scale, stage] <br>
where N is the name of the parameter,<br>
val is a Float64 for the parameter value,<br>
scale is a Float64 for multiplying the step size by (0.0 for frozen),<br>
and lastly stage is an integer.

The stage is mostly ignored by the program but you can use negative values to enforce a shared value for the operators.
For example for two equivalent tops, set the stage of F_2, ρz_2, ρx_2, and Vn_2 to -4 and each of their values to 1.0.
This wil cause them to use the same values as used for top 1 and maintain this during fitting.

$A$, $B$, and $C$ are the rotational constants corresponding to the $z$, $x$, and $y$ molecule fixed axes respectively. 
$D_{ab}$ is the off-diagonal rotational constant on the $x$ and $z$ axes for the Rho Axis Method. 
When it is included, $A$ and $B$ refer to their RAM meanings not their PAM meanings. These are all in MHz

$F$ is the internal rotor constant. It is the reduced value as is used in RAM36. 
$\rho_z$ is the coupling term between the methyl rotation and the $z$-axis angular momentum. 
It engages in the Hamiltonian as $-2\rho FP_{\alpha}N_{z}$ and $(A-F\rho^{2})N_{z}^{2}$. 
$\rho_x$ is analogous to $\rho_z$ though for the $x$-axis. 
For RAM fits, leave this as zero but it can be included for PAM fits.
The $Vn$ is the first term in the potential expansion as determined by NFOLD.

The $\epsilon$ terms are the spin-rotation coupling terms referring to axes in their subscripts. 
The term $\epsilon_{zx}$ or $\epsilon_{xz}$ is taken as the average these terms and is used for $T^2_{\pm1}(\epsilon)$. The term $T^1_{\pm1}(\epsilon)$ is no longer in the second order Hamiltonian but is available as a special higher order term.
These are transformed into spherical tensor notation upon start of the code.
These can also serve as nuclear spin-rotation coupling as they are mathematically equivalent.

The $\chi$ terms are the quadrupole terms, again referring to the axes in their subscript.
This can be used as either the electronic spin-spin in molecules of $S\ge1$ or as the nuclear electric quadrupole moment for molecules with $I\ge1$.
They are also transformed into spherical tensor form upon initializaiton.

Lastly are the spin-torsion coupling terms $\eta_z$ and $\eta_x$ which refer to the operators $S_zP_{\alpha}$ and $S_xP_{\alpha}$.
Definitions of the operators are provided in the paper.

**Caution**. The user facing second order terms are not directly used nor fit by the program.
The table below shows how the scale values in the 2nd order section are actually treated.
| User Facing   | Internal      | User Facing   | Internal      |
|:-------------:|:------------: |:-------------:|:------------: |
| A  | BK $N_z^2$           | ϵzz | T⁰₀(ϵ) $T^0_0(N,S)$ |
| B  | BN $N^2$             | ϵxx | T²₀(ϵ) $T^2_0(N,S)$ |
| C  | B± $(N_+^2 + N_-^2)$ | ϵyy | T²₂(ϵ) $T^2_{\pm2}(N,S)$ |
| ρz | rz $P_{\alpha}N_z$   | ϵzx | T²₁(ϵ) $T^2_{\pm1}(N,S)$ |
| ρx | rx $P_{\alpha}N_x$   | ϵxz | T¹₁(ϵ) $T^1_{\pm1}(N,S)$ |

### Higher order operators
These are manual coded in operators that are implemented as the anti-commutator of what the user codes in.
Each line must be in the structure of: <br>
N = [op, val, scale, stage] <br>
where N is the name of the parameter,<br>
op is a string defining the operator,<br>
val is a Float64 for the parameter value,<br>
scale is a Float64 for multiplying the step size by (0.0 for frozen),<br>
and lastly stage is an integer.
 

The stages are defined as:<br>
2 for one-top terms<br>
1 for top-top terms<br>
0 for rotation dependent terms.<br>
The code runs from high stage to low stage.
If the number of stages is higher than the stage of the operator, it will be moved into the highest allowed stage (e.i. in 2 stage mode, a stage 2 operator will be used in stage 1).
Negative values of stage can be used to add mutliple operators together and the value on the lines with the negative values will be multiplied by the value in the operator that has the non-negative stage.

I encourage you to maintain your own standard form for the other operators.
There are no checks of symmetry like in RAM36 so go wild but the user is also expected to know what they are doing. 
ere are some unicode characters for easier name: Δ, δ, Φ, ϕ, ϵ, χ, ρ, η, μ.

The current list of available operators and their meanings are given below:
E: Identy matrix \
Nz^n: $N_z^n$ \
N2^n: $N^{2n}$ \
Np^n: $N_+^n$ \
Nm^n: $N_-^n$ \
Npm^n: $N_+^n + N_-^n$

You can add your own operators by editing the baseops.jl file.
The operator functions must take the arguments of ψ<:Psi, p :: Int, and q::Int and must return a Sparse Matrix.
Use RPsi for the wavefunction type for rotational operators and TPsi for torsional.
VPsi will be added later to support vibrational operators.

### Intensity
The intensity section is used to define the molecular electric n-pole operators and their Fourier expansions with respect to a torsional coordinate.
The lines for each operator follow the form of: <br>
N = [op, val]<br>
where N is a string containing the name of the term,<br>
op is a string defining the operator, <br>
and val is the value of the parameter.



### Line File Format
The experimental transitions should be included in the molnam.lne file.
This uses a comma delimited format in the following structure:
```
Ju, Nu, Kau, Kcu, mu, Jl, Nl, Kal, Kcl, ml, freq, unc
```
$J$ is the total angular momentum (called $F$ for cases with nuclear spin) and $N$ is the angular momentum of the molecular frame.
$K_a$ and $K_c$ are the traditional assymetric top state labels and $m$ is the free rotor quantum number.
The u and l labels denote the upper and lower states respectively.
Lastly, freq is the experimental frequency and unc is the uncertainty in the line position which is used for weighting in the Levenberg-Marquadt routine.
Make sure that Ju, Jl, freq, and unc are all floating point numbers (include a decimal point) and all the rest are integers.

## Installation

These installation instructions are to the best of our knowledge, but both developers run the same Linux distribution. Please send us an email if it's not working! We need to know if we want to fix it.

### WINDOWS

There are two possible ways to install westerfit on Windows. The first one is running it natively in Windows, and the second is using the Windows Subsystem for Linux. The latter version is more robust, but the first one is likely easier.

#### Natively in Windows
1. Install [Julia](https://apps.microsoft.com/detail/9njnww8pvkmn?ocid=webpdpshare).
2. Open Julia. Press `]` to enter package mode. Enter `update` and then `add Westerfit`. The download may take some time.
3. Before using westerfit, you will need to create an input file, "molnam.inp", and a line list, "molnam.lne". The "molnam" string should be a helpful file name. Make sure these are in the same directory as each other. 
4. To use Westerfit, open your Command Line. Navigate to the directory your input "molnam.inp" is in. Enter `julia`. You should see the Julia startup.
```
using Westerfit
westerfit("molnam")
```
That should run westerfit!

#### Windows Subsystem for Linux
1. Set up [Windows Subsystem for Linux](https://www.howtogeek.com/744328/how-to-install-the-windows-subsystem-for-linux-on-windows-11/) with the Ubuntu App.
2. Open up Ubuntu and do the set-up process, including setting a username and password (you won't be able to see the password when you type it).
3. From here, move to the LINUX instructions.

Useful Windows Notes:
1. You can make a folder with the command `mkdir` followed by the name of the folder: `mkdir molecule`. Enter the folder with cd and see what's inside it with ls. You can view text files with vim, and exit vim with Esc+`:q`.
2. You can access your normal File Explorer by going to your home directory (`cd ~`) and typing `explorer.exe`.
3. You may run into issues where your files have the wrong type of linebreaks. To fix this, install dos2unix by typing `sudo apt install dos2unix` and typing in your password. Then, type `dos2unix FILENAME`, with your file name inserted. This should fix the problem.


### LINUX

1. Install [Julia](https://julialang.org/). It is recommended that you do so through [juliaup](https://github.com/JuliaLang/juliaup).
2. Run the command `julia` to enter a REPL session. Enter `]` to enter package mode. Enter `add Westerfit`.
3. Navigate to a directory in your PATH and create a file named `westerfit` containing the 3 lines below with X replaced by the number of threads you want to run on (more is better, you can also just remove the -tX altogether). You can use `which julia` to determine your exact path to Julia.
```
#!/PATH/TO/julia -tX
using Westerfit
westerfit(ARGS[1])
```
4. Run `chmod +x westerfit` to turn `westerfit` into an executeable. 
5. Try running westerfit! Enter `westerfit molnam` and it will run on molnam.inp. Make sure you're in the directory that molnam.inp is in. "molnam" can be replaced with any molecule name.

Enjoy! And feel free to reach out if anything goes wrong or you have any questions!


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

