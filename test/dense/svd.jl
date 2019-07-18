using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv
import LinearAlgebra.LAPACK: gebrd!

@testset "test svd $T by gesvd $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    matrixToTest = rand(T, size, size)

    # use the default Linear Algebra lib to calculate the right answer for testing
    right_answer = svd(matrixToTest).S
    S = right_answer

    # to test the GPU interface, one should convert the matrix data to cuda
    if interface == "GPU"
        matrixToTest = cu(matrixToTest)
    end

    # define the job U and job VT, which are required by the MAGMA lib
    jobu = 'A'
    jobvt = 'A'

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()

    # call the basic (overloaded) wrapper gesvd! for gesvd subroutines
    result = gesvd!(jobu,jobvt,matrixToTest)

    # in the result, the wanted answer lies in the second position
    s = result[2]

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    @test S ≈ s

    # else we print the detailed error info
    if !(S ≈ s)
        println("Unfortunately, the test failed.")
        println("Here is some possibly useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
    end
end


@testset "test svd $T by gesdd $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    matrixToTest = rand(T, size, size)

    # use the default Linear Algebra lib to calculate the right answer for testing
    right_answer = svd(matrixToTest).S
    S = right_answer

    # to test the GPU interface, one should convert the matrix data to cuda
    if interface == "GPU"
        matrixToTest = cu(matrixToTest)
    end

    # define the job MAGMA, which are required by the MAGMA lib
    job_magma = 'A'

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()

    # call the basic (overloaded) wrapper gesdd! for gesvd subroutines
    result = gesdd!(job_magma,matrixToTest)

    # in the result, the wanted answer lies in the second position
    s = result[2]

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    @test S ≈ s

    # else we print the detailed error info
    if !(S ≈ s)
        println("Unfortunately, the test failed.")
        println("Here is some possibly useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
    end
end

@testset "test svd $T by gebrd $interface" for T in [Float32, Float64, ComplexF32, ComplexF64], interface in ["CPU", "GPU"]

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    matrixToTest = rand(T, size, size)
    matrixToTest_copy = copy(matrixToTest)

    # use the default Linear Algebra lib to calculate the right answer for testing
    right_answer = gebrd!(matrixToTest)

    # to test the GPU interface, one should convert the matrix data to cuda
    if interface == "GPU"
        matrixToTest_copy = cu(matrixToTest_copy)
    end

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()

    # call the basic (overloaded) wrapper gesdd! for gesvd subroutines
    result = magma_gebrd!(matrixToTest_copy)

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()

    # if S is approximately equal to s, we defined then it's alright
    for index in 1:length(result)
        @test result[index] ≈ right_answer[index]
    end

    # else we print the detailed error info
    if !(1==1)
        println("Unfortunately, the test failed.")
        println("Here is some possibly useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
    end
end
