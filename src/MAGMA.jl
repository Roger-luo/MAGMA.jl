module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays
using CEnum

using LinearAlgebra: triu, tril, dot, checksquare

export magma_gels!, magma_gesvd!, magma_gesdd!, magmaInit, magmaFinalize, magma_gebrd!, libmagma, magmafunc_gpu

export magma_gesv!

# export some wrappers in clang auto-generation
export magma_init, magma_finalize

include("clang/libmagma_common.jl")
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

# include the files of subroutines
include("dense/dense.jl")

end  # modul MAGMA
