
This file exist to remind Wes of how the hell her program is supposed to work

Torsional basis set determining values:
-mcalc
-vtcalc
-vtmax

1 stage mode:
tsr matrix will be: (2J+1)(2S+1)(2mcalc+1)^length(nfold)

2 stage mode:
top-top matrix will be: (2mcalc+1)^length(nfold)
this will spit out vtcalc+1 eigenvalues
tsr matrix will be: (2J+1)(2S+1)(vtcalc+1)

3 stage mode:
one-top matrix will be: (2mcalc+1)
this will spit out vtcalc+1 eigenvalues per top
top-top matrix will be: (vtcalc+1)^length(nfold)
this will spit out vtcalc+1 eigenvalues total
tsr matrix will be: (2J+1)(2S+1)(vtcalc+1)

TORSIONAL LOGIC:

tvs is the one top wave functions
sandwich-tt check & process will be handled by final enact function

1top-two stage will follow the pattern of ntop-two stage

case 0: one top, one stage, final stage
--nothing
--length(nf)==1 && isnothing(wvs.top)

case 1: one top, two stage, first stage
--nothing
--length(nf)==1 && isnothing(tvs)

case 2: one top, two stage, second stage
--sandwich-tt
--length(nf)==1 && isnothing(tvs)

ONE TOP ALWAYS does nothing

case 3: ntop, one stage, only stage
--torsetter
--length(nf)>1 && isnothing(tvs) 

case 4: ntop, two stage, first stage 
--torsetter
--length(nf)>1 && isnothing(tvs) 

case 5: ntop, two stage, second stage
--torsetter, sandwich-tt 
--length(nf)>1 && isnothing(tvs)

LENGTH(NF)>1 && ISNOTHING(TVS) ALWAYS torsets

case 6: ntop, three stage, first stage
--nothing
--length(nf)>1 && !isnothing(tvs) && iszero(tvs)

case 7: ntop, three stage, second stage
--sandwich, torsetter
--length(nf)>1 && !isnothing(tvs) && !iszero(tvs)

case 8: ntop, three stage, final stage
--sandwich, torsetter, sandwich-tt
--length(nf)>1 && !isnothing(tvs) && !iszero(tvs)

if isone(length(psi.nf)) #cases 0,1,2
   #nothing
elseif length(psi.nf)>1 && isnothing(tvs) # cases 3,4,5
   torsetter!()
elseif length(nf)>1 && !isnothing(tvs) && iszero(tvs) # case 6
   #nothing
elseif length(nf)>1 && !isnothing(tvs) && !iszero(tvs) # case 7,8
   sandwich()
   torsetter!()
else
@warn "unexpected condition in evalulation of tor op"
end
