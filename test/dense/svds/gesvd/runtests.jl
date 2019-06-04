using Test, Random, LinearAlgebra

@testset "gesvd" begin
    include("magmaSgesvd.jl")
    include("magmaDgesvd.jl")
end
