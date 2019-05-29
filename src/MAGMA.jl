module MAGMA

export gesvd!#,magmaInit,magmaFinalize

# low-level wrappers of the MAGMA library
include("common.jl")
include("Utilities/init.jl")
include("Dense Linear Algebra/Singular Value Decomposition/gesvd/gesvd.jl")

end  # modul MAGMA
