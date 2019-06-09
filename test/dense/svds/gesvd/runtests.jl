using Test, Random, LinearAlgebra

using MAGMA: magmaInit, magmaFinalize, gesvd!

@testset "magmaSgesvd" begin
    include("magmaSgesvd.jl")
end

@testset "magmaDgesvd" begin
    include("magmaDgesvd.jl")
end
