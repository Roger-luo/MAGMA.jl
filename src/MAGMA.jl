module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays
using CEnum

# some useful functions from LinearAlgebra for some intermediate processes
using LinearAlgebra: triu, tril, dot, checksquare, chkstride1

# some useful const definitions or utility functions
export magmaTypeTuple, magmaTypeDict, magmaTypeList
export chkfinite

# export wrappers in svds
export magma_gels!, magma_gesvd!, magma_gesdd!, magmaInit, magmaFinalize, magma_gebrd!

# export wrappers in linearsystemsolver
export magma_gesv!, magma_getri!, magma_getrs!, magma_getrf!, magma_gerbt!, magma_gesv_rbt!, magma_geqrsv!
export magma_posv!, magma_hesv!, magma_sysv!

# export wrappers in factorization
export magma_geqrf!, magma_geqlf!, magma_gelqf!

# export wrappers in eigenvalues
export magma_geev!, magma_gehrd!

# export some wrappers in clang auto-generation
export magma_init, magma_finalize

#* some extra equipments
require_one_based_indexing(A...) = !has_offset_axes(A...) || throw(ArgumentError("offset arrays are not supported but got an array with index other than 1"))

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

# ! for parts of subroutines which can not do the query
macro magmafunc_nb(function_name)
	return Expr(:quote, Symbol("magma_get_", function_name, "_nb"))
end

macro magmafunc_generic(function_name, str)
	return Expr(:quote, Symbol(function_name, str))
end

const magmaTypeTuple= (Float32, Float64, ComplexF32, ComplexF64)
const magmaTypeList = ["Float32", "Float64", "ComplexF32", "ComplexF64"]
const magmaTypeList_real = ["Float32", "Float64"]
const magmaTypeList_complex = ["ComplexF32", "ComplexF64"]
const magmaTypeDict = Dict(Float32=>"s", Float64=>"d", ComplexF32=>"c", ComplexF64=>"z",
"Float32"=>"s", "Float64"=>"d", "ComplexF32"=>"c", "ComplexF64"=>"z",)

function chkfinite(A::AbstractMatrix)
    for a in A
        if !isfinite(a)
            throw(ArgumentError("matrix contains Infs or NaNs"))
        end
    end
    return true
end

# include the files of subroutines
include("dense/dense.jl")

end  # modul MAGMA
