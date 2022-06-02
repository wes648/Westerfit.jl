
using LinearAlgebra
using Symbolics

@variables A B C Dab

u1 = (1/sqrt(2)) .* [-1 0 1; 0 sqrt(2) 0; 1 0 1]
u2 = (1/sqrt(2)) .* [-1 0 0 0 1; 0 -1 0 1 0; 0 0 sqrt(2) 0 0; 0 1 0 1 0; 1 0 0 0 1]
h1p = [A+(B+C)/2 0 (B-C)/2; 0 B+C 0; (B-C)/2 0 A+(B+C)/2]
h1r = [A+(B+C)/2 -Dab/sqrt(2) (B-C)/2; -Dab/sqrt(2) B+C Dab/sqrt(2); (B-C)/2 Dab/sqrt(2) A+(B+C)/2]
h2p = [4*A+B+C 0 sqrt(6)*(B-C)/2 0 0; 0 A+(5/2)*(B+C) 0 (3/2)*(B-C) 0; sqrt(6)*(B-C)/2 0 3*(B+C) 0 sqrt(6)*(B-C)/2; 0 (3/2)*(B-C) 0 A+(5/2)*(B+C) 0; 0 0 (sqrt(6)/2)*(B-C) 0 4*A+B+C]
h2r = [4*A+B+C -3*Dab sqrt(6)*(B-C)/2 0 0; -3*Dab A+(5/2)*(B+C) -sqrt(6)*Dab/2 (3/2)*(B-C) 0; sqrt(6)*(B-C)/2 -sqrt(6)*Dab/2 3*(B+C) sqrt(6)*Dab/2 sqrt(6)*(B-C)/2; 0 (3/2)*(B-C) sqrt(6)*Dab/2 A+(5/2)*(B+C) 3*Dab; 0 0 (sqrt(6)/2)*(B-C) 3*Dab 4*A+B+C]


println("u1*h1p*u1")
println(u1*h1p*u1)
println()

println("u1*h1r*u1")
println(u1*h1r*u1)
println()

println("u2*h2p*u2")
println(u2*h2p*u2)
println()

println("u2*h2r*u2")
println(u2*h2r*u2)
println()
