module MAGMA

export gesvd!#,magmaInit,magmaFinalize

include("common.jl")
include("Utilities/init.jl")
include("Dense Linear Algebra/Singular Value Decomposition/gesvd/gesvd.jl")

end  # modul MAGMA
