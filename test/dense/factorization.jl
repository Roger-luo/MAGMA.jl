import LinearAlgebra.LAPACK

@testset "geqrf" begin
    @testset for elty in magmaTypeTuple
        @testset for interface in (Array, CuArray)
            n = 2
            m = 2
            A = interface(rand(elty, n, m))
            Atest = copy(A)

            A, tau = LAPACK.geqrf!(Matrix(A))
            Atest, tau_test = magma_geqrf!(Atest)
            @test A ≈ Atest
            @test tau ≈ tau_test
        end
    end
end

# @testset "gerqf" begin
#     @testset for elty in magmaTypeTuple
#         @testset for interface in (Array, CuArray)
#             n = 2
#             m = 2
#             A = interface(rand(elty, n, m))
#             Atest = copy(A)

#             A, tau = LAPACK.gerqf!(Matrix(A))
#             Atest, tau_test = magma_gerqf!(Atest)
#             @test A ≈ Atest
#             @test tau ≈ tau_test
#         end
#     end
# end

@testset "geqlf" begin
    @testset for elty in magmaTypeTuple
        @testset for interface in (Array, CuArray)
            n = 2
            m = 2
            A = interface(rand(elty, n, m))
            Atest = copy(A)

            A, tau = LAPACK.geqlf!(Matrix(A))
            Atest, tau_test = magma_geqlf!(Atest)
            @test A ≈ Atest
            @test tau ≈ tau_test
        end
    end
end

@testset "gelqf" begin
    @testset for elty in magmaTypeTuple
        @testset for interface in (Array, CuArray)
            n = 2
            m = 2
            A = interface(rand(elty, n, m))
            Atest = copy(A)

            A, tau = LAPACK.gelqf!(Matrix(A))
            Atest, tau_test = magma_gelqf!(Atest)
            @test A ≈ Atest
            @test tau ≈ tau_test
        end
    end
end
