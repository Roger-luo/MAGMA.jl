module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays
using CEnum

using LinearAlgebra: triu, tril, dot, checksquare

# export wrappers
export magma_gels!, magma_gesvd!, magma_gesdd!, magmaInit, magmaFinalize, magma_gebrd!

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

const magmaTypeList = ["Float32", "Float64", "ComplexF32", "ComplexF64"]
const magmaTypeDict = Dict(Float32=>"s", Float64=>"d", ComplexF32=>"c", ComplexF64=>"z",
"Float32"=>"s", "Float64"=>"d", "ComplexF32"=>"c", "ComplexF64"=>"z",)

# include the files of subroutines
include("dense/dense.jl")

end  # modul MAGMA
