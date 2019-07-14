using MAGMA
using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

const error_threshold = 1e-6

@testset "test svd $T by gesvd" for T in [Float32, Float64, ComplexF32, ComplexF64]
    matrixToTest = rand(T, 2, 2)

    right_answer = svd(matrixToTest).S
    S = right_answer

    jobu = 'A'
    jobvt = 'A'

    success=magmaInit()

    result = gesvd!(jobu,jobvt,matrixToTest)

    s = result[2]

    magmaFinalize()

    diff = S .- s
    error_value = norm(diff)

    @test error_value < error_threshold

    if error_value >= error_threshold
        println("Unfortunately, the test failed.")
        println("Here is some maybe useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
    end
end

@testset "test svd $T by gesdd" for T in [Float32, Float64, ComplexF32, ComplexF64]
    matrixToTest = rand(T, 2, 2)

    right_answer = svd(matrixToTest).S
    S = right_answer

    job_magma = 'A'

    ldu=2
    ldvt=2
    lwork=400
    success=magmaInit()

    result = gesdd!(job_magma,matrixToTest)

    s = result[2]

    magmaFinalize()

    diff = S .- s
    error_value = norm(diff)

    @test error_value < error_threshold

    if error_value >= error_threshold
        println("Unfortunately, the test failed.")
        println("Here is some maybe useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
    end
end


@testset "test svd $T by gesdd" for T in [Float32, Float64, ComplexF32, ComplexF64]
    randomGeneratedMatrix = rand(T, 2, 2)
    matrixToTest = cu(randomGeneratedMatrix)

    right_answer = svd(randomGeneratedMatrix).S
    S = right_answer

    job_magma = 'A'

    ldu=2
    ldvt=2
    lwork=400
    success=magmaInit()

    result = gesdd!(job_magma,matrixToTest)

    s = result[2]

    magmaFinalize()

    diff = S .- s
    error_value = norm(diff)

    @test error_value < error_threshold

    if error_value >= error_threshold
        println("Unfortunately, the test failed.")
        println("Here is some maybe useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
    end
end
