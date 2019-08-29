# wrappers for MAGMA factorization functions

const function_list = ("geqrf", "gerqf", "geqlf", "gelqf")
for type in magmaTypeList
    # create the symbols for element types
    elty = Symbol(type)
    # generate the symbol variables for our wrappers
    for func_name in function_list
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end
    
    @eval begin

        function magma_geqrf!(A::)

    end

end
