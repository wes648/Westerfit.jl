# westerfit
A new program for the simulating and fitting of molecular rotational spectra for $C_s$ molecules with one internal rotor and one spin source.

The pre-print of the paper is available [here](https://dx.doi.org/10.2139/ssrn.4807560).

Please feel free to direct any questions about the program to wes.harper2@gmail.com

This is me making an edit
# Quickstart

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

**νmin** (0.0::Float64): Lowest frequency in the simulation in GHz. The first character is a *nu* not a v so be careful

**νmax** (40.::Float64): Highest frequency in the simulation in GHz. Make sure this is larger than νmin or the code will crash

**INTthres** (0.0::Float64): Minimum intensity value included in the file

**TK** (8.0::Float64): Temperature of the simulation, not sure this is actually used because a certain author finds thermodynamics abhorrent

**λlm0** (0.0001::Float64): Scale-factor used to determine the inital Levenberg-Marquardt Parameter. This gets multiplied by a function of the rms to determine the LBMQ Parameter used in a given step.

**turducken** (1::Int): Number of Levenberg-Marquardt steps calculated on single step before recalculated the Jacobian. Doesn't seem worth it given the current performance of the Jacobian but can give an occasional performance boost

**assign** (expect::String): Determines how the quantum numbers are assigned after diagonalization. The default, **expect**, uses the expectation values of $m$ & then $N$ followed by an energetic sorting to assign $K$. **RAM36** uses the contributions of different blocks of the eigenvectors to provide a spin-expanded version of the assigner in RAM36. **expectk** is similar to expect but uses the expectation values of $K$ as well and doesn't seem to work. Lastly, **eeo** does the expection values of $m$ and $N$ followed by eigenvector analysis to assign $K$. This does the best job of reproducing SPFIT's DIAG=3. Generally **expect** is recommended but **RAM36** works very nicely for single perturbations (spin or torsion). I'm personally fond of the theory in **eeo** but find its performance lacking.


### Second Order Parameters
The second order Hamiltonian of westerfit is hardcoded. You can comment out any lines with # but the number of lines must remain fixed in this section. All the parameters default to zero and there are some internal transformations that occur upon initializing the code.

$A$, $B$, and $C$ are the rotational constants corresponding to the $z$, $x$, and $y$ molecule fixed axes respectively. 
$D_{ab}$ is the off-diagonal rotational constant for the Rho Axis Method. 
When it is included, $A$ and $B$ refer to their RAM meanings not their PAM meanings. These are all in MHz

$F$ is the internal rotor constant. It is the reduced value as is used in RAM36. 
$\rho_z$ is the coupling term between the methyl rotation and the $z$-axis angular momentum. 
It engages in the Hamiltonian as $-2\rho FP_{\alpha}N_{z}$ and $(A-F\rho^{2})N_{z}^{2}$. 
$\rho_x$ is analogous to $\rho_z$ though for the $x$-axis. 
For RAM fits, leave this as zero but it can be included for PAM fits.
The $Vn$ is the first term in the potential expansion as determined by NFOLD.

The $\epsilon$ terms are the spin-rotation coupling terms referring to axes in their subscripts. 
If only one of either $\epsilon_{zx}$ or $\epsilon_{xz}$ is provided, the code will internally treat the value as average term, $T^2_{\pm}(\epsilon)$.
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
| C  | B⨦ $(N_+^2 + N_-^2)$ | ϵyy | T²₂(ϵ) $T^2_{\pm2}(N,S)$ |
| ρz | rz $P_{\alpha}N_z$   | ϵzx | T²₁(ϵ) $T^2_{\pm1}(N,S)$ |
| ρx | rx $P_{\alpha}N_x$   | ϵxz | T¹₁(ϵ) $T^1_{\pm1}(N,S)$ |

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

###Line File Format
The experimental transitions should be included in the molnam.lne file.
This uses a comma delimited format in the following structure:
```
Ju, Nu, Kau, Kcu, mu, Jl, Nl, Kal, Kcl, ml, freq, unc
```
$J$ is the total angular momentum ($F$ for cases with nuclear spin) and $N$ is the angular momentum of the molecular frame.
$K_a$ and $K_c$ are the traditional assymetric top state labels and $m$ is the free rotor quantum number.
The u and l labels denote the upper and lower states respectively.
Lastly, freq is the experimental frequency and unc is the uncertainty in the line position which is used for weighting in the Levenberg-Marquadt routine.
Make sure that Ju, Jl, freq, and unc are all floating point numbers (include a decimal point) and all the rest are integers.

## Installation
Unfortuantely, westerfit installation is a bit of a mess. 
Part of this is from the dev team both running the same Linux distribution and not having other machines for testing. 
I'm afraid the WIGXJPFjl installer seems to be broken again. I'm working on it but manual installation of WIGXJPF may be necessary at the current time

westerfit is now available in the Julia package manager!

### WINDOWS

1. Set up [Windows Subsystem for Linux](https://www.howtogeek.com/744328/how-to-install-the-windows-subsystem-for-linux-on-windows-11/) with the Ubuntu App.
2. Open up Ubuntu and do the set-up process, including setting a username and password (you won't be able to see the password when you type it).
3. From here, move to the LINUX instructions.

Useful Windows Notes:
1. You can make a folder with the command `mkdir` followed by the name of the folder: `mkdir molecule`. Enter the folder with cd and see what's inside it with ls. You can view text files with vim, and exit vim with Esc+`:q`.
2. You can access your normal File Explorer by going to your home directory (`cd ~`) and typing `explorer.exe`.
3. You may run into issues where your files have the wrong type of linebreaks. To fix this, install dos2unix by typing `sudo add dos2unix` and typing in your password. Then, type `dos2unix FILENAME`, with your file name inserted. This should fix the problem.


### LINUX

1. Install [Julia](https://julialang.org/). It is recommended that you do so through [juliaup](https://github.com/JuliaLang/juliaup).
2. Make sure you have a C compiler, such as gcc. You can check this by typing `which gcc` into your terminal.
3. Run the command `julia` to enter a REPL session. Enter `]` to enter package mode. Enter `add Westerfit`.
4. Unfortunately, WIGXJPFjl, a dependency, is causing issues. You can try to remedy this by typing `add WIGXJPFjl`. If that doesn't work, you can try manual compilation and putting the libwigxjpf_shared.so library in the package's deps directory; please let me know if you do this. If you don't feel comfortable trying that, feel free to shoot me an email.
5. Open your bashrc with the command `vim ~/.bashrc`. Press `i` and enter to enter edit mode. Add these three lines, with X replaced by the number of threads you want to run on (more is better, you can also just remove the -tX altogether):
```
#!/PATH/TO/julia -tX
using Westerfit
westerfit(ARGS[1])
```
Then enter `. ~/.bashrc` to update your terminal session.
6. Try running westerfit! Enter `westerfit molnam` and it will run on molnam.inp. Make sure you're in the directory that molnam.inp is in.

To install westerfit, you will need to install [Julia](https://julialang.org/) and I recommend doing so through [juliaup](https://github.com/JuliaLang/juliaup).
You also need to have a C compiler available for the WIGXJPF dependency.
Ideally, once Julia is installed, type `julia` into your command line to enter a REPL session and hit `]` to enter package mode.
Simply type `add Westerfit` and you should be good to go!

The automated installation will add the SparseArrays, StaticArrays, and WIGXJPFjl packages to your Julia installation.
The WIGXJPFjl package is a wrapper for [WIGXJPF](http://fy.chalmers.se/subatom/wigxjpf/) and is currently the main point of failure in the install process.
It should be smoothed up thanks to the BinaryBuilder.jl based install script but I've only tested it on Linux.
You may have to manually compile and place the libwigxjpf_shared.so library in the package's deps directory and I ask that you let me know if that's the case.
You may also have better luck using the Ubuntu terminal from the Microsoft Store on Windows if you don't have access to a Linux machine.
The code also uses DelimitedFiles, LinearAlgebra, Printf, and Dates but I believe all of those are included in Base.

If the automated installation method fails, please shoot me an email.
I'll gladly help you troubleshoot the installation and hopefully fix it so it installs smoothly for the next person.

Lastly, I recommend making a simple runner script called `westerfit` somewhere in your path.
This script just needs the following few lines:
```
#!/PATH/TO/julia -tX
using Westerfit
westerfit(ARGS[1])
```
Set X to be the number of threads you want westerfit to run on (more is better!) and fill in your full path to the julia.
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

