# sample code for magma binding
export MagmaAllVec, MagmaIVec, MagmaNoVec


# MAGMA enum constants
# the whole file will be stored as enums.jl
#just like in JuliaGPU/MAGMA.jl

# MAGMA constants indicating the vectors status
# as input/output for some functions
# For example, the gesvd functions will use
# MagmaNoVec, MagmaSomeVec, MagmaAllVec and
# MagmaOverwriteVec to indicate the
# strategies that will be applied to the SVD
# U matrix and VT matrix in A = U Σ V**T
# (for MagmaOverwriteVec it is going to overwrite A)
const MagmaNoVec         = 301
const MagmaVec           = 302
const MagmaIVec          = 303
const MagmaAllVec        = 304
const MagmaSomeVec       = 305
const MagmaOverwriteVec  = 306
const MagmaBacktransVec  = 307

### some certain setting for this testing
using CUDAdrv
using CUDAapi
using CUDAnative
using CuArrays
# indicate the library position
const libmagma = "/usr/local/magma/lib/libmagma.so"
