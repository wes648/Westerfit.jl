
using WIGXJPFjl,LinearAlgebra
Δlist(J::Real,S::Real) = (abs(J-S)):(J+S)
powneg1(k::Real)::Int = isodd(k) ? -1 : 1
δi(x::Real,y::Real)::Int = x==y
□rt(x::Real)::Float64 = √(x*(x>zero(x)))

#atom A
a = 1.5
χzz = -27665
χxx =   -418
χyy =  28083
χxz =  49499
χa = [χzz; -√(2/3)*χxz; √(1/6)*(χxx-χyy)]
#atom b
b = 1.0
χzz =   1297
χxx =  -2796
χyy =   1499
χxz =  -1088
χb = [χzz; -√(2/3)*χxz; √(1/6)*(χxx-χyy)]


qna = [fill(a,Int((2a+1)*(2b+1))) kron(collect(-a:a),ones(Int(2b+1)))]
qnb = [fill(b,Int((2a+1)*(2b+1))) kron(ones(Int(2a+1)),collect(-b:b))]
qns = reduce(vcat, [[fill(s,Int(2s+1)) collect(-s:s)] for s ∈ Δlist(a,b)])

coupling(s,Σ,a,α,b,β) = powneg1(a-b+Σ)*√(2s+1)*wig3j(a,b,s, α,β,-Σ)
coupling(s,a,b) = coupling.(s[:,1],s[:,2], a[:,1]',a[:,2]', b[:,1]',b[:,2]')

function hq(χ1,χ2,χ3,a,αb,βb,αk,βk)
   out  = 0.5*χ1*(3*αk^2 - a*(a+1))*δi(αb,αk)*δi(βb,βk)
   out += χ2*(αk+0.5)*□rt(a*(a+1) - αk*(αk+1))*δi(αb,αk+1)*δi(βb,βk)
   out += χ2*(αk-0.5)*□rt(a*(a+1) - αk*(αk-1))*δi(αb,αk-1)*δi(βb,βk)
   out += χ3*□rt(a*(a+1) - αk*(αk+1))*□rt(a*(a+1) - (αk+1)*(αk+2))*δi(αb,αk+2)*δi(βb,βk)
   out += χ3*□rt(a*(a+1) - αk*(αk-1))*□rt(a*(a+1) - (αk-1)*(αk-2))*δi(αb,αk-2)*δi(βb,βk)
   return out
end
hq(χ,a,α,β) = hq.(χ[1],χ[2],χ[3],a,α',β',α,β)

W = coupling(qns,qna,qnb)
H = hq(χa,a,qna[:,2],qnb[:,2]) + hq(χb,b,qnb[:,2],qna[:,2])