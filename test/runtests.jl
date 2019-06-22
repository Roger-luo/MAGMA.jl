using MAGMA
using Test, Random, LinearAlgebra

@testset "dense linear algebra" begin
    include("dense/svd.jl")
end
