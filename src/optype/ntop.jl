
function σcount(nfold::Real)::Int
   if isodd(nfold)
      out = ceil(Int,0.5*nfold)
   elseif iseven(nfold)&&nfold≠0
      out = floor(Int,nfold/4) + 1
   else #nfold == 0
      out = 1
   end
   return out
end
function msgen(T::Type,nfold::Real,mcalc::Real,σ::Real)
   if nfold==0
      return ones(T,1)
   else
   lim = mcalc*nfold
   if (σ==0)&&(isodd(nfold))
      marray = collect(T,-lim:nfold:lim)
   elseif (σ==0)&&(iseven(nfold))
      lim = floor(Int,lim/2)
      marray = collect(T,-lim:floor(T,nfold/2):lim)
   elseif (σ≠0)&&(iseven(nfold))
      lim = floor(Int,lim/2)
      marray = collect(T,-lim+σ:floor(T,nfold/2):lim+σ)
   else #(σ≠0)&&(isodd(nfold))
      marray = collect(T,(-lim+σ):nfold:(lim+σ))
   end
   return marray
   end
end
function msgen(nfold::Int,mcalc::Int,σ::Int)::Array{Int}
   marray = msgen(Int,nfold,mcalc,abs(σ))
   if σ < 0
      marray .*= -1
   end
   return marray
end

function σgen(nf1::Int,nf2::Int)::Array{Int,2}
   σcnt1 = σcount(nf1)
   σcnt2 = σcount(nf2)-1
	out = zeros(Int,2,0)
   for i in 1:σcnt1
		σ1 = i-1
		if σ1==0
			σ2 = collect(0:σcnt2)'
		else
			σ2 = collect(σcnt2:-1:-σcnt2)'
		end#if
		out = hcat(out, [fill(σ1,length(σ2))'; σ2])
	end#for
	return out
end#function

function σ2gen(σ1,nf2)
   σcnt2 = σcount(nf2)-1
   if σ1==0
      σ2 = collect(0:σcnt2)
   else
      σ2 = collect(σcnt2:-1:-σcnt)
   end
   out = hcat(fill(σ1,length(σ2)), σ2)
   return out
end
function nextσ0(nfc,nf2::Int)::Array{Int,2}
   σcnt2 = σcount(nf2)-1
   σ2 = collect(0:σcnt2)
   t = zeros(1,nfc)
   out = hcat(repeat(t,length(σ2)), σ2)
   return out'
end
function nextσ1(σs::Array,nf2::Int)::Array{Int,2}
   σcnt2 = σcount(nf2)-1
   σ2 = collect(σcnt2:-1:-σcnt2)
   out = hcat(repeat(σs',length(σ2)), σ2)
   return out
end
function σgen_indef(nf::Array{Int})::Array{Int,2}
   old = σgen(nf[1],nf[2])
   for i in 3:length(nf)
      σlsti = nf[i] - iseven(nf[i])
      new = vcat(old,fill(0,size(old,2))')
      for j in 2:σlsti
         σ = Int(cospi(j)*floor(j/2))
         part = vcat(old[:,1+(j==σlsti):end],fill(σ,size(old[:,1+(j==σlsti):end],2))')
         new = hcat(new,part)
      end#for j
      old = new
   end#for i
   return old
end
function σgen_indef(nf::Int)::Array{Int}
   return collect(0:σcount(nf[1])-1)'
end


function msgen_2top(nf::Array,mcalc)::Array{Int,3}
   σs = σgen(nf[1],nf[2])
   out = zeros(Int,2mcalc+1,length(nf),size(σs,1))
   for i in 1:size(σs,2), j in 1:length(nf)
      out[:,j,i] = msgen(nf[j],mcalc,σs[i,j])
   end
   return out
end

function msgen_all_indef(nf::Array,mcalc::Int)::Array{Int,3}
   σs = σgen_indef(nf)
   out = zeros(Int,2mcalc+1,length(nf),size(σs,1))
   for i in 1:size(σs,1), j in 1:length(nf)
      out[:,j,i] = msgen(nf[j],mcalc,σs[i,j])
   end
   return out
end
function msgen_indef(nf::Array{Int},mcalc::Int,σs::Array{Int})::Array{Int,2}
   out = zeros(Int,2mcalc+1,length(nf))
   for j in 1:length(nf)
      out[:,j] = sort!(msgen(nf[j],mcalc,σs[j]))
   end
   return out
end

function ntopset!(nf,ms,tid,mat::SparseMatrixCSC{Float64,Int64})::SparseMatrixCSC{Float64,Int64}
   lm = size(ms,1)
   ln = length(nf)
   t = lm*(tid - 1)
   if t≠0
      mat = kron(I(lm*(tid - 1)), mat)
   end
   t = lm*(ln-tid)
   if t≠0
      mat = kron(mat, I(lm*(ln-tid)) )
   end
   return mat
end
function ntopset!(nf,ms,tid,mat::Array)::Array
   lm = size(ms,1)
   ln = length(nf)
   t = lm*(tid - 1)
   if t≠0
      mat = kron(ones(lm*(tid - 1)), mat)
   end
   t = lm*(ln-tid)
   if t≠0
      mat = kron(mat, ones(lm*(ln-tid)) )
   end
   return mat
end

