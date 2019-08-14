using MAGMA
using Test, Random, LinearAlgebra

using Test, LinearAlgebra, CuArrays, CUDAnative, CUDAapi, CUDAdrv

@testset "dense linear algebra" begin
    include("dense/svd.jl")
    include("dense/linearsystemsolver.jl")
end
