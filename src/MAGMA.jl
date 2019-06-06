module MAGMA

export gesvd!, libmagma#,magmaInit,magmaFinalize

include("common.jl")
include("utils/Utilities.jl")
include("dense/svd/gesvd.jl")

end  # modul MAGMA
