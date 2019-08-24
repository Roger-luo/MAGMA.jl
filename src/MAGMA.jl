module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays
using CEnum

using LinearAlgebra: triu, tril, dot, checksquare

export magma_gels!, magma_gesvd!, magma_gesdd!, magmaInit, magmaFinalize, magma_gebrd!, libmagma, magmafunc_gpu

export magma_gesv!

# export some wrappers in clang auto-generation
export magma_init

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
# include("enums.jl")
include("clang/libmagma_common_v2.jl")
include("clang/libmagma_api.jl")



"""
the path to magma binary
"""
const libmagma = "/usr/local/magma/lib/libmagma.so"

macro magmafunc(function_name)
	return Expr(:quote, Symbol("magma_", function_name))
end

macro magmafunc_gpu(function_name)
	return Expr(:quote, Symbol("magma_", function_name, "_gpu"))
end

macro magmatype(elty)
	if elty == Float32
		return 's'
	end
	if elty == Float64
		return 'd'
	end
	if elty == ComplexF32
		return 'c'
	end
	if elty == ComplexF64
		return 'z'
	end
end

# >>> The following are some Utility functions' wrappers >>>
# magma_init
# function magmaInit()
# 	success = ccall((:magma_init, libmagma),Cint,())
# 	if success != 0
# 		println("MAGMA initiation error with success = ", success)
# 	end
# end

# magma_finalize
function magmaFinalize()
	ccall((:magma_finalize, libmagma),Cint,())
end

# <<< End of wrappers for Utility

# include the files of subroutines
include("dense/dense.jl")

end  # modul MAGMA
