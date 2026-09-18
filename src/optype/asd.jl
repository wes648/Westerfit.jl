"""
This is the Automated Suggestor of Distortions

user input:
[asd]
t1 = ["A", a]
t2 = ["B", b]
t3 = ["C", c]

the code will then generate expectation value vectors of:
Δ⟨AA⟩
Δ⟨AB⟩
Δ⟨AC⟩
Δ⟨BB⟩
Δ⟨BC⟩
Δ⟨CC⟩

function std(vec::Vector)::Float64
   mean = sum(vec)/length(vec)
   √(sum(abs2,vec - mean)/length(vec))
end

they will then be normalized and compared to the normalized omc vector
tAA = std(omc/N - Δ⟨AA⟩/M)
...

if std is small operator is suggested
if mean = 0
   param_suggest: M/N

"""