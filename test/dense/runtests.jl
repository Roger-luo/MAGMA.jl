using MAGMA
using Test, LinearAlgebra, Random

@testset "singular_value_decomposition" begin
    include("svds/runtests.jl")
end
