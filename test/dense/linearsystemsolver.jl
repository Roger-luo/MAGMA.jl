using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv
import LinearAlgebra.LAPACK: gels!, gesv!

@testset "test magma_gels! $T by gels $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    A = rand(T, size, size)
    B = rand(T, size, size)

    dA = cu(A)
    dB = cu(B)

    A_copy = copy(A)
    B_copy = copy(B)

    (A_test, B_test) = interface == "GPU" ? (dA, dB) : (A_copy, B_copy)

    right_answer = gels!('N', A, B)

    magmaInit()
    # result = interface == "GPU" ? magma_gels!(MAGMA.MagmaNoTrans, dA, dB) : magma_gels!(MAGMA.MagmaNoTrans, A_copy, B_copy)
    result = magma_gels!(MAGMA.MagmaNoTrans, A_test, B_test)
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    for i in 1:length(result)
        @test result[i] ≈ right_answer[i]
    end

end

@testset "test magma_gesv! $T by gesv $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size_test = 2
    A = rand(T, size_test, size_test)
    B = rand(T, size_test, size_test)

    A = one(A)
    B = one(B)

    dA = cu(A)
    dB = cu(B)

    A_copy = copy(A)
    B_copy = copy(B)

    (A_test, B_test) = interface == "GPU" ? (dA, dB) : (A_copy, B_copy)

    right_answer = gesv!(A, B)
    println("Right answer is: \n", right_answer)

    magmaInit()
    # result = magma_gesv!(A_test, B_test)
    n = LinearAlgebra.checksquare(A)
    isGPU = (A isa CuArray && B isa CuArray)
    # if size(B,1) != n
    #     throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
    # end
    lda  = max(1,stride(A,2))
    ldb  = max(1,stride(B,2))
    ipiv = similar(Matrix(A), Int32, n)
    info = Ref{Cint}()
    ccall((:magma_sgesv, libmagma), Cint,
          (Ref{Cint}, Ref{Cint},
           PtrOrCuPtr{Float32}, Cint, Ptr{Cint},
           PtrOrCuPtr{Float32}, Cint, Ptr{Cint}),
           n, size(B,2),
           A, lda, ipiv,
           B, ldb, info)
    magmaFinalize()
    result = (B,A,ipiv)



    # if S is approximately equal to s, we defined then it's alright
    for i in 1:length(result)
        @test result[i] ≈ right_answer[i]
    end

end
