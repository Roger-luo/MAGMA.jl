using MAGMA
using Test, Random, LinearAlgebra

# for simply tests the very beginning codes for gesvd
@testset "MAGMA.jl" begin
    # Write your own tests here.
    include("Dense Linear Algebra/runtests.jl")
end
