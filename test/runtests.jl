using MAGMA
using Test, Random, LinearAlgebra

@testset "dense linear algebra" begin
    include("dense/svd.jl")
    # include("dense/linearsystemsolver.jl")
end

@testset "Debug gels" begin
    using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv
    import LinearAlgebra.LAPACK: gels!

    T = Float32

    # randomly generate a 2 by 2 matrix for testing
    size = 2
    A = rand(T, size, size)
    B = rand(T, size, size)

    dA = cu(A)
    dB = cu(B)

    # initialize the MAGMA lib, serving as a necessary part before working
    magmaInit()

    println("CPU interface for gels!")
    println(magma_gels!(0, A, B))

    magma_gels!(0, dA, dB)

    # finalize the MAGMA lib, serving as a necessary part after working
    magmaFinalize()
end
