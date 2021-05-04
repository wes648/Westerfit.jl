
module WIGXJPF

function __init__()
    max_two_j = 1000
    ccall((:wig_table_init, "./libwigxjpf_shared.so"),
          Cvoid,
          (Cint, Cint),
          max_two_j, 9)
    ccall((:wig_temp_init, "./libwigxjpf_shared.so"),
          Cvoid,
          (Cint,),
          max_two_j)
end

doubled(i::Integer) = 2i
doubled(f::Float64) = Int(2f)
doubled(r::Rational) = Int(2r)

function wig3jj(j12::Integer, j22::Integer, j32::Integer,
                m12::Integer, m22::Integer, m32::Integer)
    ccall((:wig3jj, "./libwigxjpf_shared.so"),
          Cdouble,
          (Cint, Cint, Cint,
           Cint, Cint, Cint),
          j12, j22, j32,
          m12, m22, m32)
end

wig3j(j1, j2, j3,
      m1, m2, m3) = wig3jj(doubled(j1), doubled(j2), doubled(j3),
                           doubled(m1), doubled(m2), doubled(m3))

function wig6jj(j12::Integer, j22::Integer, j32::Integer,
                j42::Integer, j52::Integer, j62::Integer)
    ccall((:wig6jj, "./libwigxjpf_shared.so"),
          Cdouble,
          (Cint, Cint, Cint,
           Cint, Cint, Cint),
          j12, j22, j32,
          j42, j52, j62)
end

wig6j(j1, j2, j3,
      j4, j5, j6) = wig6jj(doubled(j1), doubled(j2), doubled(j3),
                           doubled(j4), doubled(j5), doubled(j6))

function wig9j(j12::Integer, j22::Integer, j32::Integer,
               j42::Integer, j52::Integer, j62::Integer,
               j72::Integer, j82::Integer, j92::Integer)
    ccall((:wig9jj, "./libwigxjpf_shared.so"),
          Cdouble,
          (Cint, Cint, Cint,
           Cint, Cint, Cint,
           Cint, Cint, Cint),
          j12, j22, j32,
          j42, j52, j62,
          j72, j82, j92)
end

wig9j(j1, j2, j3,
      j4, j5, j6,
      j7, j8, j9) = wig9jj(doubled(j1), doubled(j2), doubled(j3),
                           doubled(j4), doubled(j5), doubled(j6),
                           doubled(j7), doubled(j8), doubled(j9))

export wig3jj, wig3j, wig6jj, wig6j, wig9jj, wig9j

end # module
