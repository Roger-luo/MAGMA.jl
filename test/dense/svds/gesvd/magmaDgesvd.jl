using MAGMA

using Test, LinearAlgebra

@testset "random double precision matrices" begin

    matrixToTest = rand(Cdouble, 2, 2)

    right_answer = svd(matrixToTest).S
    S = right_answer

    jobu = 'A'
    jobvt = 'A'

    ldu=2
    ldvt=2
    lwork=400
    success=magmaInit()

    result = gesvd!(jobu,jobvt,matrixToTest)

    s = result[2]

    magmaFinalize()

    diff = S .- s
    error_value = norm(diff)

    @test error_value < 1e-7

    if error_value >= 1e-7
        println("Unfortunately, the test failed.")
        println("Here is some maybe useful information:")
        println("the element_type is ", typeof(matrixToTest))
        println("the right answer = ", S)
        println("However, MAGMA got the answer = ", s)
        println("The info returned by MAGMA is: ", info[1])
    end

end
