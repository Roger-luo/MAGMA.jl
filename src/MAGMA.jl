module MAGMA

export gesvd!, libmagma#,magmaInit,magmaFinalize

include("common.jl")
include("Utilities/Utilities.jl")

include("Dense Linear Algebra/Singular Value Decomposition/gesvd/gesvd.jl")

end  # modul MAGMA
