# westerfit
A new program for the simulating and fitting of molecular rotational spectra for C_s molecules with one internal rotor and one spin source.

Description of current status:
Torsion-Rotation tested against BELGI & RAM36. See test/belgi and example/2ba
Hyper-fine tested against SPFTI. see test/spfit
SR undergoing testing. see test/spfit

The westerfit input file is divided into three sections each with a designated header. The very first line is a title card and the sections are the Control Settings (%CNTRLS), Second Order Parameters (%2NDORDER), and the additional Operators & Parameters (%PARAMS).

## Control Settings
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


