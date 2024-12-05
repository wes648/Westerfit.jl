
function opdictinit()::Dict{SubString,Op}
dict = Dict{String,Op}("Nz" => Nz, "N2" => N2, "Np" => Np, "Nm" => Nm, 
   "Npm" => Npm,"Nx" => Nx, "sNy" => sNy, "NS" => NS, "S2" => S2, "Sz" => Sz, 
   "Sp" => Sp, "Sm" => Sm, "Spm" => Spm, "Sx" => Sx, "sSy" => sSy, "Pα" => Pα, 
   "cosα" => cosα, "Pβ" => Pβ, "cosβ" => cosβ, "Pγ" => Pγ, "cosγ" => cosγ)
   return dict
end

function snip(s)
   if s[1]==' '
      return s[2:end]
   else
      return s
   end
end
