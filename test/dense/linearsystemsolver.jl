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

    magma_init()
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

    dA = cu(A)
    dB = cu(B)

    A_copy = copy(A)
    B_copy = copy(B)

    (A_test, B_test) = interface == "GPU" ? (dA, dB) : (A_copy, B_copy)

    right_answer = gesv!(A, B)

    magma_init()
    result = magma_gesv!(A_test, B_test)
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    for i in 1:length(result)
        @test result[i] ≈ right_answer[i]
    end

end
