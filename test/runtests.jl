using MAGMA
using Test, Random, LinearAlgebra

# for simply tests the very beginning codes for gesvd
@testset "dense_linear_algebra" begin
    # Write your own tests here.
    include("dense/runtests.jl")
end
