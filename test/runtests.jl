using MAGMA, LinearAlgebra.LAPACK
using Test, Random, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

@testset "dense linear algebra" begin
    include("dense/svd.jl")
    include("dense/linearsystemsolver.jl")
    include("dense/factorization.jl")
    include("dense/eigenvalues.jl")
end
