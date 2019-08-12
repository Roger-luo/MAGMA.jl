using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv
import LinearAlgebra.LAPACK: gels!

@testset "test gels $T by gesvd $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    A = rand(T, size, size)
    B = rand(T, size, size)

    dA = cu(A)
    dB = cu(B)

    A_copy = copy(A)
    B_copy = copy(B)
    #
    # println("Before gels")
    # println("A = ", A)
    # println("B = ", B)

    # use the default Linear Algebra lib to calculate the right answer for testing
    right_answer = gels!('N', A, B)
    #
    # println("After gels")
    # println("A = ", A)
    # println("B = ", B)

    # to test the GPU interface, one should convert the matrix data to cuda
    # if interface == "GPU"
    #     matrixToTest = cu(matrixToTest)
    # end


    # call the basic (overloaded) wrapper gesvd! for gesvd subroutines
    if interface == "GPU"
        # initialize the MAGMA lib, serving as a necessary part before working
        magmaInit()
        result = magma_gels!(MAGMA.MagmaNoTrans, dA, dB)

        # finalize the MAGMA lib, serving as a necessary part after working
        magmaFinalize()
    else
        # initialize the MAGMA lib, serving as a necessary part before working
        magmaInit()

        result = magma_gels!(MAGMA.MagmaNoTrans, A_copy, B_copy)

        # finalize the MAGMA lib, serving as a necessary part after working
        magmaFinalize()
    end


    # if S is approximately equal to s, we defined then it's alright
    for i in 1:length(result)
        @test result[i] ≈ right_answer[i]
    end

end
