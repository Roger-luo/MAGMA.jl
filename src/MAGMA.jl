module MAGMA
using CUDAdrv, CUDAapi, CUDAnative, CuArrays

export gesvd!, gesdd!, magmaInit, magmaFinalize

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
       set_vector, set_matrix, get,
       axpy!, dot,
       gemv!

import Base: length,size


# some types' definitions and their accordingly chararcters
const magmaTypes = Union{Float32, Float64, ComplexF32, ComplexF64}
const magmaTypesList =  (Float32,Float64,ComplexF32,ComplexF64)
const typeChar = Dict(
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
struct set_vector_struct{T <: magmaTypes}
    ptr::Ptr{T}
    length::Int
end

struct set_matrix_struct{T <: magmaTypes}
    ptr::Ptr{T}
    size::Tuple{Int, Int}
end

length(x::set_vector_struct) = x.length
size(x::set_matrix_struct)   = x.size
size(x::set_matrix_struct,i) = x.size[i]


include("enums.jl")

for T in magmaTypesList

	# communication between CPU and GPU
    setvector = ("magma_$(typeChar[T])setvector", libmagma)
    getvector = ("magma_$(typeChar[T])getvector", libmagma)
    setmatrix = ("magma_$(typeChar[T])setmatrix", libmagma)
    getmatrix = ("magma_$(typeChar[T])getmatrix", libmagma)

    axpy = ("magma_$(typeChar[T])axpy",libmagma)
    dot  = T <: Complex ? ("magma_$(typeChar[T])dotc",libmagma) :
                         ("magma_$(typeChar[T])dot", libmagma)
    gemv = ("magma_$(typeChar[T])gemv",libmagma)

    @eval function set_vector(x::Vector{$T})
        ptrptr = Array(Ptr{$T},1)
        err = ccall((:magma_malloc, libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Errorlibmagmadevice memory ($err)")
        ptr = ptrptr[1]
        ccall($setvector,Cvoid,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x,1,ptr,1,"","",0)
        set_vector_struct{$T}(ptr,length(x))
    end

    @eval function get_vector(x::set_vector_struct{$T})
        y = Array($T,length(x))
        ccall($getvector,Cvoid,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x.ptr,1,y,1,"","",0)
        y
    end

    @eval function set_matrix(x::Matrix{$T})
        ptrptr = Array(Ptr{Float64},1)
        err = ccall((:magma_malloc,libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Error allocating device memory ($err)")
        ptr = ptrptr[1]
        ccall($setmatrix,Cvoid,(Cint,Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               size(x,1),size(x,2),x,size(x,1),ptr,size(x,1),"","",0)
        set_matrix{$T}(ptr,size(x))
    end

    @eval function get_matrix(x::set_matrix_struct{$T})
        y = Array($T,length(x))
        ccall($getmatrix,Cvoid,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x.ptr,1,y,1,"","",0)
        y
    end


    # Level 1 BLAS

    @eval function axpy!(a::$T,x::set_vector_struct{$T},y::set_vector_struct{$T})
        length(x) != length(y) && error("Vectors must be the same length")
        ccall($axpy,Cvoid,(Cint,$T,Ptr{$T},Cint,Ptr{$T},Cint),
                          length(x),a,x.ptr,1,y.ptr,1)
        nothing
    end

    @eval function dot(x::set_vector_struct{$T},y::set_vector_struct{$T})
        length(x) != length(y) && error("Vectors must be the same length")
        ccall($dot,$T,(Cint,Ptr{$T},Cint,Ptr{$T},Cint),
                       length(x),x.ptr,1,y.ptr,1)
    end

    # Level 2 BLAS

    @eval function gemv!(a::$T,A::set_matrix_struct{$T},x::set_vector_struct{$T},b::$T,y::set_vector_struct{$T})
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

# include the files of subroutines
include("dense/dense.jl")

end  # modul MAGMA
