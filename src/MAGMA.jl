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

# https://github.com/mweastwood/MAGMA.jl

export device_reset,
       MagmaVector, MagmaMatrix, get,
       axpy!, dot,
       gemv!

import Base: length,size,dot


# some types' definitions and their accordingly chararcters
const magmatypes = Union{Float32, Float64, ComplexF32, ComplexF64}
const magmatypeslist =  (Float32,Float64,ComplexF32,ComplexF64)
const typechar = Dict(
    Float32 => 's',
    Float64 => 'd',
    ComplexF32 => 'c',
    ComplexF64 => 'z'
)

"""
the path to magma binary
"""
const libmagma = "/usr/local/magma/lib/libmagma.so"


# https://github.com/mweastwood/MAGMA.jl
struct MagmaVector{T <: magmatypes}
    ptr::Ptr{T}
    length::Int
end

struct MagmaMatrix{T <: magmatypes}
    ptr::Ptr{T}
    size::Tuple{Int, Int}
end

length(x::MagmaVector) = x.length
size(x::MagmaMatrix)   = x.size
size(x::MagmaMatrix,i) = x.size[i]


include("enums.jl")

for T in magmatypeslist
    setvector = ("magma_$(typechar[T])setvector_internal", libmagma)
    getvector = ("magma_$(typechar[T])getvector_internal", libmagma)
    setmatrix = ("magma_$(typechar[T])setmatrix_internal", libmagma)
    axpy = ("magma_$(typechar[T])axpy",libmagma)
    dot  = T <: Complex ? ("magma_$(typechar[T])dotc",libmagma) :
                         ("magma_$(typechar[T])dot", libmagma)
    gemv = ("magma_$(typechar[T])gemv",libmagma)

    @eval function MagmaVector(x::Vector{$T})
        ptrptr = Array(Ptr{$T},1)
        err = ccall((:magma_malloc, libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Errorlibmagmadevice memory ($err)")
        ptr = ptrptr[1]
        ccall($setvector,Cvoid,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x,1,ptr,1,"","",0)
        MagmaVector{$T}(ptr,length(x))
    end

    @eval function get(x::MagmaVector{$T})
        y = Array($T,length(x))
        ccall($getvector,Cvoid,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x.ptr,1,y,1,"","",0)
        y
    end

    @eval function MagmaMatrix(x::Matrix{$T})
        ptrptr = Array(Ptr{Float64},1)
        err = ccall((:magma_malloc,libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Error allocating device memory ($err)")
        ptr = ptrptr[1]
        ccall($setmatrix,Cvoid,(Cint,Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               size(x,1),size(x,2),x,size(x,1),ptr,size(x,1),"","",0)
        MagmaMatrix{$T}(ptr,size(x))
    end

    # Level 1 BLAS

    @eval function axpy!(a::$T,x::MagmaVector{$T},y::MagmaVector{$T})
        length(x) != length(y) && error("Vectors must be the same length")
        ccall($axpy,Cvoid,(Cint,$T,Ptr{$T},Cint,Ptr{$T},Cint),
                          length(x),a,x.ptr,1,y.ptr,1)
        nothing
    end

    @eval function dot(x::MagmaVector{$T},y::MagmaVector{$T})
        length(x) != length(y) && error("Vectors must be the same length")
        ccall($dot,$T,(Cint,Ptr{$T},Cint,Ptr{$T},Cint),
                       length(x),x.ptr,1,y.ptr,1)
    end

    # Level 2 BLAS

    @eval function gemv!(a::$T,A::MagmaMatrix{$T},x::MagmaVector{$T},b::$T,y::MagmaVector{$T})
        size(A,2) != length(x) && error("Matrix and vectors have inconsistent sizes")
        size(A,1) != length(y) && error("Matrix and vectors have inconsistent sizes")
        ccall($gemv,Cvoid,(Cint,Cint,Cint,$T,Ptr{$T},Cint,Ptr{$T},Cint,$T,Ptr{$T},Cint),
                          MagmaNoTrans,size(A,1),size(A,2),a,A.ptr,size(A,1),x.ptr,1,b,y.ptr,1)
        nothing
    end
end
# end of fork https://github.com/mweastwood/MAGMA.jl

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
function magmaMalloc_CPU()
	ccall()
end
# <<< End of wrappers for Utility

include("dense/dense.jl")


end  # modul MAGMA
