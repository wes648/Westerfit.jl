tsrdiag -- torsion spin rotation diagonalizer
initializes wang transformation
checks whether sigma is zero
	sigma is the symmetry flag for the torsion -- if it is divisible by three (A state) it is zero, if it is not divisible by three (E state) it is + or - 1 (here we only cover +1 but they're degenerate).
	If sigma is zero, we can do a second layer of the Wang transformation on the torsions
ur is a Wang transformation of the rotational part, ut is a wang transformation of the torsional part

H..mat function builds the Hamiltonian matrix
H, in either case, builds the matrix and then performs the Wang transformation. it also makes ti a dense matrix instead of a sparse matrix

ln 72
jacobisweep fxn does a fixed number of iterations of the jacobi diagonalization approach -- takes two args (matrix and interation count). returns jacobi pre-conditioned (part of the diagonalization) (H) and rvec (rotational vectors, product of all rotational matrices)
LAPACK.syev! diagonalizes the pre-conditioned matrix -- faster than jacobi but scrambles the eigenvectors. 
assign goes throug the vectors Lapack made and finds the dominant element in each column, arranges them into an identity-like matrix. then combines with the rvecs to get the true eigenvectors. it also generates the quantum numbers. returns qn, eigenvalues, eigenvectors

tsrdiag is the central fxn of westerfit and westersim -- they have different fxns wrapped around it, but this is where the eigenvalues and vectors come from

TSRCALC -- torsion spin rotation (energy level) calculator
first thing it does it determine the minimum possible J value -- if S is half integer, it's 1/2, and if S is a whole number, it's 0
jmax is found from nmax - s (will be fixed)
jarray -- builds an array of all the j values
jd -- calculates the j degeneracies. 
outvals, outvecs, outqns -- creates empty arrays for the eigenvalues, eigenvectors, and quantum numbers. may be lies.

@threads is the parallelized calling of tsrdiag -- goes through every value in jarray
	finds start index and final index using jinds (j indices)
	tqns, etc (temporary) -- calls tsrdiag for every j
	si, fi -- finds the limit of mmax in the temporary stack of quantum numbers
	then we put our tvals into outvals in the spots where they go as calculated by sind/find and si/fi, etc with outvecs and outqns

WESTERSIM -- the simulator

sigmacount -- derives number of sigmas from NFOLD
	same initializations as TSRCALC except uses sigma count as an additional dimenion on the arrays (will be removed from tsrcalc because seriously why)

line 119 -- calculates qns, vals, vecs for sigma = 0 only, writes to first sets of the things
transitions -- calculates the transitions (slow function)
then we iterate through all the other sigmas (b/c transitions needs to already exist to concatenate w/ temptrns)

JLISTER -- finds unique J and sigma pairs in the input line list
LIMEIGCALC -- limited eigen calculator. only calculates tsrdiag for j-sigma pairs
DOGLEG -- only being kept for --memes-- "a potentially better method of optimization"

WESTERFIT_HANDCODED -- the fitter, doesn't work for NFOLD>1
currently defines mmax, mcalc, S for test parameters, defines tsrparams and lines instead of calling
lineprep -- determines indicies for the lines, grabs the frequencies and uncertainties. input for lines: quantum ns, frequency, uncertainty
use jlister
define nmax
scales -- not fully implemented, binary float/notfloat
lambda -- levenberg-marquardt parameter (placeholder)
lbmq_opttr -- leverberg-marquardt optimizer trust region, runs entire LM fitting routine (has a bad name because it does not use a trust region actually, but it could be turned back on or separately callable so who knows)
output difference between fit results and original parameters

WESTERENERGIES (not in use rn)

pred2lne -- takes westersim transition output and converts it to westerfit line input format (b/c there's often extra information in a simulation)

perturbation randomizes your input a little based on scale and floating

FILE: FILEHANDLING

LINEPREP -- converts the input line list into a code firendly format, already documented
EngWriter -- energy writer, produces two files of all the different energy levels, one for A states one for E states. convrets the energies into wavenumbers and then builds and prints strings to the file, comma delimited, looks fixed width

TraWriteSPCAT -- transition writer spcat, write an spcat version of the file so that it can interact with pickett stuff (such as using ABES), already well documented
TraWriter -- produces a comma delimited value structure, makes the outputs look prettier

File: JACOBI, well-documented by Christina already, should be adapted to sparse arrays to increase speed

File: ASSIGN

ASSIGN -- assigns the quantum numbers of the eigenvalues
	run leadfact to create assigner (the permutation)
	then perform assigner on the values
	then perform assigner on the syev vectors and multiply by the rotation vectors from jacobi
	if statement for diabatic sorter (not used rn)
	QNs -- qngen runs
	returns QNs, vals, vecs, all in proper orders
	

LEADFACT -- leading factor, done on a copy of vecs to prevent it from becoming zeroes. finds the permutation to change the order post-syev into the correct order.
	turns the vectors into absolute values (only looking at magnitude).
	build the array equla to number of states in length filled with zeroes (c, named thus for unknown reasons)
	for each state, find the largest element in the vecs, returns the index (tuple b/c it's a square matrix)
	pull out the state number -- which of the basis states it is
	pull out eigenvalue number -- eigenvalue from energetic sorting
	zero out that column and row so we don't double assign it
	sorts the eigenvalue numbers to the states (kth lowest energy level is associated with teh basis step with index s)
	sortperm -- generates the permutation to put the eigenvalues in the correct order

Adiabatic sorter -- use the same state labeling even when the energies flip places, not sure it works rn

finds every state of the same symmetry, abs(k), abs(m) and then energetically sort them with their labels attached

other shit that is not in use

File: Common. already well-documented and mostly understandable anyway.

on 2n+1 for degeneracy -- "it's not accurate, but it stems from a place of accuracy"

"if this reads as a little bit deranged, it's because it is" ~on paramprep

File: OPTIMIZER

RMSCALC -- RMS calculator. takes the calculated values and the indicies of the transitions and the observed frequencies (ofreqs)
	cfreqs -- calculated frequencies generation from energies
	omc is the difference between calculated and actual
	rms is an rms duh

ANADERIV -- analytical derivative
	takes the derivative with respect to a given paramaeter by setting the others to zero and that one to one
	Follows hellmann-feynman theorem
	do a wang transformation on the derivative hamiltonian matrix (b/c built with rp)
	use the eigenvectors to calculate the expectation value

LBMQ_GAIN -- not currently being used, method of assessing whether to expand or shrink a trust region

build_jcbn! -- jcbn gets passed in, rewritten, then sent out
	iterates over every index and finds their matched vectors 
	iterates over the permutation (list of floated parameters) -- don't make a derivative for a parameter we're not floating

build_hess! -- makes the hessian. also makes the transpose of the jacobian times the weights

approx2dirdrv! -- not being used because it is a fit of lunacy, also includes lbmq_acc!

lbmq_step! -- simplest implementation of LBMQ
	make a hermissian matrix of the hessian plus the parameter times the diagonal of the hessian (gets rid of all off-diag elements by setting to zero). core to the LBMQ method. Lambda (parameter) controls the ratio of the two things -- when it gets small, we behave like gauss-newton, when it's big, we behave like steepest descent. Use diagonal Hessian instead of identity to maintain scaling of the parameters.
	while isposdef -- forces A to be a positive definite matrix -- every eigenvalue would be positive and unique and you can do a Cholesky decomposition on it, which allows for ldiv to run faster (left division). start w/ smaller lambda, then double it until we get a pos def matrix. provides a large amount of stability b/c cholesky can fail when too different initial guess
	cholesky! -- matrix decomposition. does stuff. necessary and improves performance of ldiv
	ldiv -- left division -- finds the step size by comparing A to the negative gradient
	returns beta (step -- array of values by which parameters will be shifted) and lambda (LBMQ parameter)

lbmq_turducken! -- three LBMQ steps, all a little bit different, also because memes
	first repeat lmbq_step
	want to keep track of individual steps and also total sum of step -- hence why betafull exists
	shift the floated parameters
	recalculate vals and vectors
	recalculate rms based on new values -- to get omc
	use new OMC to calculate a new step
	repeat w/o calculating the jacobian (duck)
	do it again w/o calculating the jacobian again (chicken)
three steps off of one calculation of the jacobian so it all runs faster (fixes some of the overperturbation issues common to lbmq)

paramunc! -- calculates parameter uncertainty (diagnonal elements of inverse of Hessian)

lbmq_opttr -- LMBQ optimizer trust region (previously discussed as being poorly named)
	calculate values and vectors
	store old versions in separate variable name
	do rmscalc
	ddetermine the permutation -- make scales into a sparse vector, then find non-zero elements
	goal -- want to work in later, add something based on statistically allowed lowest RMS (prevents articifally precise data)
	W -- finds the weights by uilding a diagonal matrix based on inverse uncertainty values
	epsilons -- convergence thresholds, one for change in RMS and the other for change in step
	LIMIT is max number of iterations
	lambdalm is LMBQ parameter (later will called from user input)
	deltalm is trust region, which gets changed by parameters floated
	counter and stoit -- not being used, but the idea is to stop it from converging on high RMS by doing a stochastic corrections
	make a whole bunch of empty matrices
	build jacobian and hessian

	while loop -- begin the lbmq process
	do the turducken step
	do a check -- find rms shift percentage, if new rms is lower, print step size, redefine new values to main name, rebuild the jacobian and Hessian, increase counter, shrink LMBQ parameter by a factor of three (if it gets too small,  ake it zero)
	if the new rms is higher, up the LMBQ parameter and repeat
	
	if RMS is less than goal, print CFOUR joke and end program #memes
	if RMS has stopped decreasing but is still higher than goal, say such has happened and end program
	if the step size if effectively zero, say such, and end program
	if the counter is over the limit, say so, and break
	some other condition, continue

	print uncertainties at the end, returns params and values

FILE: Hamiltonian

theta and phi are coupling terms originally used by Raynes. included special versions for s=0.5 to reduce floating pt operations

fh and gh are condensations of some common terms

Ediag is the energies of the diagonal matrix elements. currently not being used.

Htorm0 is the on-diagonal in m component of the torsional hamiltonian. currently fairly minimal for test purposes, but can be replaced at some point with a more complex set of parameters.
Htorm1 is the off-diagonal by 1 m element. includes a zero element b/c it would get confused about how to use an array to build a matrix where all the values are constant
Hr0K, Hr1K, Hr2K (hamiltonian rotational off-diag in K).
Hrot -- builds rotational hamiltonian. If N is zero, gives 1x1 sparse array w/ nothing in it to because it will freak out if it does not have off-diag elements and does not have long enough arrays.
	Generate the k array of all possible k values, then figure out their dimensions. (replace all degeneracy above w/ dimension).
	Generate the three arrays
	use sparse diagm to make an nd by nd matrix, assign arrays to correct locations

Htor -- torsional hamiltonian for a closed-shell molecule
	tnp -- nd, but with a special name because Wes wasn't there at the time
	if mcalc is zero, returns 0 matrix
	tmp -- md
	generate an array of k's and an array of m's -- use kron (kronicker products) to make everything the right size
	make your arrays, builds the sparse matrix

Htor -- radical version
	same as other in general
	ks is done a little differently because it has to figure out what n is from j and s
	matrix building, yadda yadda yadda

Htr -- torsion rotation hamiltonian
	makes the matrix
	rewrites the matrix as kronicker matrix of sparse identity matrix and rotational hamiltonian matrix -- gives your diagonal line of rotational hamiltonian matrix
	add the torsional matrix to it

Hsthings -- spin rotation elements following Rayne's. does not currently include the brown and sears centrifugal distortion terms (although it can at a time)

hqelem and hqmat are not being used because the WIGXJPF isn't parallelized

Hspi -- spin hamiltonian
	see above

(old) Htsr0N and 1n only implement p_α S_z operator at this time
        make karray and marray
        calculate initial array of elements
        make correct length, spdiagm
Hsr -- spin rotation Hamiltonian
        open with srprep to get all possible N and dimensions and indices and j dimension
        make the matrix, fill in the first N block
        add the three in the bottom right corner, then again the next three in the bottom right corner, etc

        Htsr -- old, not in use. [originally put the torsion on N with the m-levels inside of that, then add sr later. When adding hyperfine or a second spin, realized it made more sense to build the F blocks, have multiple F blocks for each N value (here J value)]

Htsr0Nv and 1Nv are the variants that are currently in use
"I'm suddenly feeling very bad for you"

Htsrmat2 follows new tsr structure
        has srpart and tspart, which cover the same Ns and are built in the same style similarly to sr Hamiltonian
        make an m array, do kronicker product to make things the right size, then add everything together

FILE: WIGXJPF

__init__ -- initialization of the package, sets up some functions
doubled functions -- return integer of twice the input value

wig3jj -- uses a 3j symbol where every element has been doubled
wig3j -- calls wig3jj and gives it doubled values
etc for 6j and 9j symbols

FILE: intensities, aka tracalc

T -- very sloppy, takes the zxy array of μ's and converts it to spherical tensor element q

intelem -- intensity element
        out -- checks the wigner symbol and the tensor component existence
        if it is not zero, it does the rest of the calculation -- iterates over every single element

intmat -- intensity matrix
        figures out the size of jd bra and jd ket
        mat and omat are because the dipole moments are best represented as a fourier series dependent on dihedral angle -- omat handles the off-diagonal in m terms
        make our matrices and add them together
        If sigma is zero, we do a two-layered wang transformation, if it isn't, we do a single layer -- gets better agreement between this matrix and the eigenvectors

intcalc -- not in use

partitioncalc --
        qs = partition functions. calculate that for every state, then add them all up and return them
thermal factor --
        move the partition function to a thermal factor

tracalc --
        set some things up
        run intmat
        for every combination of states
            finding Δj to make sure it is ≤2, if it is we continue
            assign jjb and jk, recalculate the lengths (because they have random zeroes at the end)
            calculate the frequency
            if either of the j values don't match the old j values, recalculate the mu matrix
            if the frequency is greater than the minimum and less than the maximum (0-80 GHz rn, but later will be user input)
                take expectation value of intensity and square it, multiply by thermal factor
            if intensity is greater than the intensity threshold, we write a line to a big array that contains the information
            if the frequency is less than the minimum and greater than the negative of the maximum, we do similar things but instead write the -freq to the array
            if the frequency is outside of those ranges, we do nothing b/c it's outside of the range of the instrument

part of the reason this is slow is because trans does not have a size beforehand
also because it has to run umat as it goes, not doing all of the umat and then calling them in parallel from a non-wig fxn

pred2lne -- converts westersim to westerfit
everything else is not in use
