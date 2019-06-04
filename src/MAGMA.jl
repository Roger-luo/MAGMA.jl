module MAGMA

export gesvd!, libmagma#,magmaInit,magmaFinalize

include("common.jl")
include("utilities/Utilities.jl")
include("dense/Dense.jl")

include("dense/svds/gesvd/gesvd.jl")

end  # modul MAGMA
