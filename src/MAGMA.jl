module MAGMA

export device_reset,
       MagmaVector, MagmaMatrix, get,
       axpy!, dot,
       gemv!

import Base: length,size,dot

using  CUDArt
import CUDArt: device_reset

const magmatypes = Union{Float32, Float64, Complex64, Complex128}
const magmatypeslist =  (Float32,Float64,Complex64,Complex128)
const typechar = Dict(
    Float32 => 's',
    Float64 => 'd',
    Complex64 => 'c',
    Complex128 => 'z'
)

type MagmaVector{T <: magmatypes}
    ptr::Ptr{T}
    length::Int
end

type MagmaMatrix{T <: magmatypes}
    ptr::Ptr{T}
    size::Tuple{Int, Int}
end

length(x::MagmaVector) = x.length
size(x::MagmaMatrix)   = x.size
size(x::MagmaMatrix,i) = x.size[i]

include("enums.jl")
const libmagma = Libdl.find_library("libmagma", ["/usr/local/magma/lib"])

for T in magmatypeslist
    setvector = ("magma_$(typechar[T])setvector_internal", libmagma)
    getvector = ("magma_$(typechar[T])getvector_internal", libmagma)
    setmatrix = ("magma_$(typechar[T])setmatrix_internal", libmagma)
    axpy = ("magma_$(typechar[T])axpy",libmagma)
    dot  = T <: Complex? ("magma_$(typechar[T])dotc",libmagma) :
                         ("magma_$(typechar[T])dot", libmagma)
    gemv = ("magma_$(typechar[T])gemv",libmagma)

    @eval function MagmaVector(x::Vector{$T})
        ptrptr = Array(Ptr{$T},1)
        err = ccall((:magma_malloc, libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Errorlibmagmadevice memory ($err)")
        ptr = ptrptr[1]
        ccall($setvector,Void,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x,1,ptr,1,"","",0)
        MagmaVector{$T}(ptr,length(x))
    end

    @eval function get(x::MagmaVector{$T})
        y = Array($T,length(x))
        ccall($getvector,Void,(Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               length(x),x.ptr,1,y,1,"","",0)
        y
    end

    @eval function MagmaMatrix(x::Matrix{$T})
        ptrptr = Array(Ptr{Float64},1)
        err = ccall((:magma_malloc,libmagma),Cint,(Ptr{Ptr{$T}},Cint),ptrptr,length(x)*sizeof($T))
        err != 0 && error("Error allocating device memory ($err)")
        ptr = ptrptr[1]
        ccall($setmatrix,Void,(Cint,Cint,Ptr{$T},Cint,Ptr{$T},Cint,
                               Ptr{UInt8},Ptr{UInt8},Cint),
                               size(x,1),size(x,2),x,size(x,1),ptr,size(x,1),"","",0)
        MagmaMatrix{$T}(ptr,size(x))
    end

    # Level 1 BLAS

    @eval function axpy!(a::$T,x::MagmaVector{$T},y::MagmaVector{$T})
        length(x) != length(y) && error("Vectors must be the same length")
        ccall($axpy,Void,(Cint,$T,Ptr{$T},Cint,Ptr{$T},Cint),
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
        ccall($gemv,Void,(Cint,Cint,Cint,$T,Ptr{$T},Cint,Ptr{$T},Cint,$T,Ptr{$T},Cint),
                          MagmaNoTrans,size(A,1),size(A,2),a,A.ptr,size(A,1),x.ptr,1,b,y.ptr,1)
        nothing
    end
end

end
