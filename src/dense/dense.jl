function char_to_magmaInt(option::AbstractChar)
    if option == 'A'
            return MagmaAllVec
    elseif option == 'S'
            return MagmaSomeVec
    elseif option == 'O'
            return MagmaOverwriteVec
    elseif option == 'N'
            return MagmaNoVec
    end
end

subsetrows(X::AbstractVector, Y::AbstractArray, k) = Y[1:k]
subsetrows(X::AbstractMatrix, Y::AbstractArray, k) = Y[1:k, :]

include("svd.jl")
include("linearsystemsolver.jl")
