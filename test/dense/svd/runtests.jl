using Test, Random, LinearAlgebra

@testset "magmaSgesvd" begin
    include("magmaSgesvd.jl")
end

@testset "magmaDgesvd" begin
    include("magmaDgesvd.jl")
end
