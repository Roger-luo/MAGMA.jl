@testset "geev" begin
    # complex is easier for now
    @testset for elty in MAGMA.magmaTypeTuple
        A = rand(elty,10,10)
        magma_result = MAGMA.magma_geev!('N','V',copy(A))
        lapack_result= LAPACK.geev!('N', 'V', copy(A))
        for i in 1:length(magma_result)
            println("Testing the $(i)th result: ")
            @test magma_result[i] ≈ lapack_result[i]
        end
    end
end
