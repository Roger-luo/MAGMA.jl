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

    result = magma_gels!(MAGMA.MagmaNoTrans, A_test, B_test)

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

    result = magma_gesv!(A_test, B_test)

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
            A = LAPACK.getri!(A, ipiv)

            iB = inv(B)
            B, ipivB = LAPACK.getrf!(B)
            B = interface(B)
            B = magma_getri!(B, ipivB)
            
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
            B = getrs!(trans_char, A, ipiv, B)
            Atest, ipivtest = magma_getrf!(Atest)
            Btest = magma_getrs!(trans, Atest, ipivtest, Btest)
            @test B ≈ Btest
        end
    end
end

@testset "getrf" begin
    @testset for elty in MAGMA.magmaTypeTuple
        @testset for interface in (Array, CuArray)
            # println(magmaTypeTuple)
            A = interface(rand(elty, 2, 2))
            B = copy(Matrix(A))

            A, ipiv = magma_getrf!(A)
            B, ipivB,infoB= LAPACK.getrf!(B)

            @test A ≈ B
            @test ipiv ≈ ipivB
        end
    end
end

@testset "gesv_rbt" begin
    @testset for elty in (Float32, Float64)#MAGMA.magmaTypeTuple
        local n = 2
        A = Array(rand(elty, n, n))
        B = Array{elty}(I, n, n)

        A_backup = Matrix(A)
        B_backup = Matrix(B)
        result = inv(Matrix(A))

        A, B, info = magma_gesv_rbt!(A, B, MAGMA.MagmaTrue)
        @test (A_backup * Matrix(B)) ≈ B_backup
    end
end

@testset "posv" begin
    @testset for elty in MAGMA.magmaTypeTuple
        @testset for interface in (CuArray, Array)
            local n = 10
            A = rand(elty,n,n)/100
            A += real(diagm(0 => n*real(rand(elty,n))))
            if elty <: Complex
                A = A + A'
            else
                A = A + transpose(A)
            end
            B = rand(elty,n,n)

            Dtest = interface(copy(A))
            Ctest = interface(copy(B))
            Dtest, Ctest, info = magma_posv!(MAGMA.MagmaUpper, Dtest, Ctest)
            
            @test A\B ≈ Ctest
        end
    end
end

@testset "hesv" begin
    @testset for elty in (ComplexF32, ComplexF64)
        # @testset for interface in (Array)
        A = rand(elty,10,10)
        A = A + A' #hermitian!
        b = rand(elty,10)
        c = A \ b
        # A = interface(A)
        # b = interface(b)
        A, b, iter, info = magma_hesv!(MAGMA.MagmaUpper, A, b)
        @test b ≈ c
        # println("Iter = $(iter)\nInfo = $(info)")

        # b,A = LAPACK.hesv!('U',A,b)
        # @test b ≈ c
        # @test_throws DimensionMismatch LAPACK.hesv!('U',A,rand(elty,11))
        # A = rand(elty,10,10)
        # A = A + A' #hermitian!
        # b = rand(elty,10)
        # c = A \ b
        # b,A = LAPACK.hesv_rook!('U',A,b)
        # @test b ≈ c
        # @test_throws DimensionMismatch LAPACK.hesv_rook!('U',A,rand(elty,11))
        # end
    end
end
