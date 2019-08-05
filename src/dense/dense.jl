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


include("svd.jl")
include("linearsystemsolver.jl")
