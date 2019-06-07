module Dense

using CUDAdrv, CuArrays, CUDAapi, CUDAnative

import MAGMA: libmagma

export gesvd!

include("svds/SVD.jl")

end
