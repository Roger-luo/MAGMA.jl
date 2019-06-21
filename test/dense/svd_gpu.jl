using MAGMA
using Test, LinearAlgebra, CuArrays

const error_threshold = 1e-6

@testset "test svd GPU $T" for T in [Float32, Float64, ComplexF32, ComplexF64]

    matrixToTest = rand(T, 2, 2)
    matrixToTest = cu(matrixToTest)

    right_answer = svd(matrixToTest).S
    S = right_answer

    jobu = 'A'
    jobvt = 'A'

    ldu=2
    ldvt=2
    lwork=400
    success=magmaInit()

    result = gesvd_gpu!(jobu,jobvt,matrixToTest)

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
