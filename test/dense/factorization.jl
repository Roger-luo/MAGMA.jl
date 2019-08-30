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
