# westerfit
A new program for the simulating and fitting of molecular rotational spectra

The full second order Hamiltonian is coded in.

The new version uses a single diagonalization stage in the free rotor quantum number for better treatment of the group theory, especially in the A-states.

The simulator works nearly in full. Statistical mechanics and spin-statistical weights need to be included but a predicted spectrum that seems believable can be produced.

The fitter is in progress. There is an option to use Optim.jl. I don't love interfacing with the Optim package so I'm working on the handcoded optimizer. Currently that does not work well but does sort of work.


Usage of the current state ofthe  program will require adjusting the paths. It is centered in tsr_fit.jl which currently generates a simulated spectrum and then perturbes the parameters and tries to fit it. That file also contains a "dehyperfined" version of the Hirota data set.
