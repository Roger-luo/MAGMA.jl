using MAGMA
using Test, Random, LinearAlgebra

using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

@testset "dense linear algebra" begin
    include("dense/svd.jl")
    # include("dense/linearsystemsolver.jl")
end

@testset "Debug gels" begin


    T = Float32

    # randomly generate a 2 by 2 matrix for testing
    A = rand(T, 2, 2)
    B = rand(T, 2, 2)

    dA = cu(A)
    dB = cu(B)


    right_answer = magma_gels!(0, A, B)
    println("CPU interface for gels! is as below:")
    println(right_answer)

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()


    m = Cint(2)
    n = Cint(2)
    info  = Ref{Cint}()
    work  = Vector{T}(undef, 1)
    lwork = Cint(-1)

    for i = 1:2  # first call returns lwork as work[1]
        # Only for debugging
        # dA = CuPtr{$elty}()
        # dB = CuPtr{$elty}()
        # println("\n","dA = ", dA, "\n")

        # print(libmagma)

        ccall((:magma_sgels3_gpu, libmagma), Cint,
              (Cint, Cint, Cint, Cint,
               CuPtr{Float32}, Cint, CuPtr{Float32}, Cint,
               Ptr{Float32}, Cint, Ptr{Cint}),
              MAGMA.MagmaNoTrans, m, n, size(B,2),
              dA, max(1,stride(A,2)), dB, max(1,stride(B,2)),
              work, lwork, info)
        # println(dA, dB)

        # chkmagmaerror(info[])
        if i == 1
            lwork = ceil(Int, real(work[1]))
            resize!(work, lwork)
        end
    end

    k   = min(m, n)
    F   = m < n ? tril(dA[1:k, 1:k]) : triu(dA[1:k, 1:k])
    ssr = Vector{T}(undef, size(dB, 2))
    for i = 1:size(dB,2)
        x = zero(T)
        for j = k+1:size(dB,1)
            x += abs2(dB[j,i])
        end
        ssr[i] = x
    end
    subsetrows(X::AbstractMatrix, Y::AbstractArray, k) = Y[1:k, :]
    result = (F, subsetrows(dB, dB, k), ssr)

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()

    for i in 1:length(right_answer)
        println("The ",i," element calculated by CPU is ", right_answer[i])
        println("The ",i," element calculated by GPU is ", result[i])
        @test right_answer[i] ≈ result[i]
    end
end
