
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
   else
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

function σgen(nf1,nf2)::Array{Int,2}
   σcnt1 = σcount(nf1)
   σcnt2 = σcount(nf2)-1
	out = zeros(Int,0,2)
   for i in 1:σcnt1
		σ1 = i-1
		if σ1==0
			σ2 = collect(0:σcnt2)
		else
			σ2 = collect(σcnt2:-1:-σcnt2)
		end#if
		out = vcat(out, [fill(σ1,length(σ2)) σ2])
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
   out = hcat(fill(σ1,legnth(σ2)) σ2)
   return out
end
function σgen2(nf::Array{Int})::Array{Int,2}
   σcnt1 = σcount(nf[1])
   out = zeros(Int,0,length(nf))
   for i in 1:σcnt1
      σ1 = i - 1
      for j in 2:length(nf)
         σ2 = σ2gen(σ1,nf[j])
      end
   end
end

function msgen(nf::Array,mcalc)::Array{Int,3}
   σs = σgen(nf[1],nf[2])
   out = zeros(Int,2mcalc+1,length(nf),size(σs,1))
   for i in 1:size(σs,1), j in 1:length(nf)
      out[:,j,i] = msgen(nf[j],mcalc,σs[i,j])
   end
   return out
end
