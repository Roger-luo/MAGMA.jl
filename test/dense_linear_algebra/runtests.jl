using MAGMA
using Test, LinearAlgebra, Random

@testset "singular_value_decomposition" begin
    include("singular_value_decomposition/runtests.jl")
end
