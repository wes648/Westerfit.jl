
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
