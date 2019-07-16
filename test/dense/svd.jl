using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

const error_threshold = 1e-6

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

    # calculate the difference between the standard answer and the calculated answer
    diff = S .- s
    error_value = norm(diff)

    # if the error value is less than the threshold we defined then it's alright
    @test error_value < error_threshold

    # else we print the detailed error info
    if error_value >= error_threshold
        println("Unfortunately, the test failed.")
        println("Here is some possibly useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
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

    # calculate the difference between the standard answer and the calculated answer
    diff = S .- s
    error_value = norm(diff)

    # if the error value is less than the threshold we defined then it's alright
    @test error_value < error_threshold

    # else we print the detailed error info
    if error_value >= error_threshold
        println("Unfortunately, the test failed.")
        println("Here is some possibly useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
    end
end
