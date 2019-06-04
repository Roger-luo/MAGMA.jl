module MAGMA

export gesvd!, libmagma#,magmaInit,magmaFinalize

include("common.jl")
include("utilities/Utilities.jl")
include("dense_linear_algebra/Dense.jl")

include("dense_linear_algebra/singular_value_decomposition/gesvd/gesvd.jl")

end  # modul MAGMA

end # module
