#!/home/wes/.julia/bin/julia

using Plots
using DelimitedFiles
using LaTeXStrings
#ys are file[:,12]
#filter by file[:,10]
#lower Ns are file[:,7]
#ofreq is file[:,10]

function lnplot(molnam::String,lns,lka,omc,rms)
   scatter(lns,omc,mz=lka,mc=cgrad(:roma,rev=true),label=false,xlabel=L"Lower $N$",
      ylabel="Residuals (kHz)",cb=:right,cbtitle=L"Lower $K_{a}$",dpi=500)
   hline!([-rms,rms],label="RMSE (kHz)",color=palette(:default)[2])
   savefig("$molnam-n.png")
end
function lnplot(molnam::String,lnsa,lkaa,omca,lnse,lkae,omce,rms)
   scatter(lnsa,omca,mz=lkaa,mc=cgrad(:roma,rev=true),label=L"$A$ States",
      xlabel=L"Lower $N$",shape=:utriangle,
      ylabel="Residuals (kHz)",cb=:right,cbtitle=L"Lower $K_{a}$",dpi=500)
   scatter!(lnse,omce,mz=lkae,label=L"$E$ States",mc=cgrad(:roma,rev=true),
      shape=:dtriangle)
   hline!([-rms,rms],label="RMSE (kHz)",color=palette(:default)[2])
   savefig("$molnam-n.png")
end

function frplot(molnam::String,fre,lka,omc,rms)
   scatter(fre,omc,mz=lka,mc=cgrad(:roma,rev=true),label=false,
       xlabel="Observed Frequency (GHz)",ylabel="Residuals (kHz)",
       cb=:right,cbtitle=L"Lower $K_{a}$",dpi=500)
   hline!([-rms,rms],label="RMSE (kHz)")
   savefig("$molnam-nu.png")
end
function frplot(molnam::String,frea,lkaa,omca,free,lkae,omce,rms)
   scatter(frea,omca,mz=lkaa,mc=cgrad(:roma,rev=true),label=L"$A$ States",
      xlabel="Observed Frequency (GHz)",shape=:utriangle,
      ylabel="Residuals (kHz)",cb=:right,cbtitle=L"Lower $K_{a}$",dpi=500)
   scatter!(free,omce,mz=lkae,label=L"$E$ States",mc=cgrad(:roma,rev=true),
      shape=:dtriangle)
   hline!([-rms,rms],label="RMSE (kHz)",color=palette(:default)[2])
   savefig("$molnam-nu.png")
end


function notor(molnam,rms)
   file = readdlm("$molnam.res",';')
   fre = file[:,11] .* 1e-3
   omc = file[:,12] .* 1e3
   lns = Int.(file[:,7])
   lka = Int.(file[:,8])
   lnplot(molnam,lns,lka,omc,rms)
   frplot(molnam,fre,lka,omc,rms)
end
function yetor(molnam,rms)
   file = readdlm("$molnam.res",';')
   fre = file[:,11] .* 1e-3
   omc = file[:,12] .* 1e3
   lns = Int.(file[:,7])
   lka = Int.(file[:,8])
   σs = Int.(file[:,10])
   frea = fre[σs.==0]
   omca = omc[σs.==0]
   lnsa = lns[σs.==0]
   lkaa = lka[σs.==0]
   free = fre[σs.==1]
   omce = omc[σs.==1]
   lnse = lns[σs.==1]
   lkae = lka[σs.==1]
   lnplot(molnam,lnsa,lkaa,omca,lnse,lkae,omce,rms)
   frplot(molnam,frea,lkaa,omca,free,lkae,omce,rms)
end

function main(molnam,torflag,rms)
   rms = parse(Float64,rms)
   if torflag=="t"
      yetor(molnam,rms)
   else
      notor(molnam,rms)
   end
end

main(ARGS[1],ARGS[2],ARGS[3])
