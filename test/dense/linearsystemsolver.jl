import LinearAlgebra.LAPACK: gels!, gesv!, getrs!, getri!

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
    magma_finalize()

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
    magma_finalize()

    # if S is approximately equal to s, we defined then it's alright
    for i in 1:length(result)
        @test result[i] ≈ right_answer[i]
    end

end

@testset "getrf/getri" begin
    @testset for elty in (Float32, Float64, ComplexF32, ComplexF64)
        @testset for interface in (CuMatrix, Matrix)
            A = rand(elty,2,2)
            B = (copy(A))

            iA = inv(A)
            A, ipiv = LAPACK.getrf!(A)
            # println("A = ", A)
            # println("ipiv = ", ipiv)
            # println("inverse of A = ", iA)
            A = LAPACK.getri!(A, ipiv)

            iB = inv(B)
            B, ipivB = LAPACK.getrf!(B)

            # println("B = ", B)
            # println("ipiv = ", ipiv)
            # println("inverse of B = ", iB)
            magma_init()
            B = interface(B)
            B = magma_getri!(B, ipivB)
            magma_finalize()
            
            @test A ≈ B
        end
    end
end

@testset "getrs" begin
    @testset for elty in (Float32, Float64, ComplexF32, ComplexF64)
        @testset for interface in (CuMatrix, Matrix)
            A = rand(elty,2,2)
            B = rand(elty,2,2)
            trans_char = 'N'

            Atest = copy(A)
            Btest = copy(B)
            trans = MAGMA.MagmaNoTrans

            A, ipiv = LAPACK.getrf!(A)
            # println("ipiv is ", ipiv)
            # println("A = ", A)
            # println("ipiv = ", ipiv)
            # println("inverse of A = ", iA)
            B = getrs!(trans_char, A, ipiv, B)

            # println("B = ", B)
            # println("ipiv = ", ipiv)
            # println("inverse of B = ", iB)
            Atest, ipivtest = LAPACK.getrf!(Atest)
            # println("ipiv of test matrix is ", ipivtest)
            magma_init()
            Btest = magma_getrs!(trans, Atest, ipivtest, Btest)
            magma_finalize()
            # println("A ≈ Atest: ", A≈Atest)
            # println("ipiv", ipiv ≈ ipivtest)
            @test_broken B ≈ Btest
        end
    end
end

@testset "getrf" begin
    @testset for elty in MAGMA.magmaTypeTuple
        @testset for interface in (Array, CuArray)
            # println(magmaTypeTuple)
            A = rand(elty, 2, 2)
            B = copy(A)

            magma_init()
            A, ipiv = magma_getrf!(A)
            magma_finalize()
            B, ipivB,infoB= LAPACK.getrf!(B)

            @test A ≈ B
            @test ipiv ≈ ipivB
        end
    end
end
