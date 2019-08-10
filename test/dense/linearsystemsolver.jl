using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv
import LinearAlgebra.LAPACK: gels!

@testset "test gels $T by gesvd $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    matrixToTest_A = rand(T, size, size)
    matrixToTest_B = rand(T, size, size)

    matrixToTest_A_copy = copy(matrixToTest_A)
    matrixToTest_B_copy = copy(matrixToTest_B)

    # use the default Linear Algebra lib to calculate the right answer for testing
    right_answer = gels!('N', matrixToTest_A, matrixToTest_B)

    # to test the GPU interface, one should convert the matrix data to cuda
    # if interface == "GPU"
    #     matrixToTest = cu(matrixToTest)
    # end

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()

    # call the basic (overloaded) wrapper gesvd! for gesvd subroutines
    result = magma_gels!(MAGMA.MagmaNoTrans, matrixToTest_A_copy, matrixToTest_B_copy)

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    @test right_answer ≈ result

end
