module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays

export gesvd!, magmaInit, magmaFinalize

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

"""
the path to magma binary
"""
const libmagma = "/usr/local/magma/lib/libmagma.so"

macro magmafunc(function_name)
    return Expr(:quote, Symbol("magma_", function_name))
end

# >>> The following are some Utility functions' wrappers >>>
# magma_init
function magmaInit()
	ccall((:magma_init, libmagma),Cint,())
end

# magma_finalize
function magmaFinalize()
	ccall((:magma_finalize, libmagma),Cint,())
end

# magma_malloc_cpu
#function magmaMalloc_CPU()

#end
# <<< End of wrappers for Utility

include("dense/dense.jl")


end  # modul MAGMA
