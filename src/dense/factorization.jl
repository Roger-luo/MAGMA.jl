# wrappers for MAGMA factorization functions

const function_list_factorization = ("geqrf", "gerqf", "geqlf", "gelqf")
for type in magmaTypeList
    # create the symbols for element types
    elty = Symbol(type)
    # generate the symbol variables for our wrappers
    for func_name in function_list_factorization
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end
    
    @eval begin

        function magma_geqrf!(A::Matrix{$elty})
            m, n  = size(A)
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)
            info  = Ref{Cint}()
            lda   = max(1, m)
            tau   = similar(A, min(m, n))

            func  = eval(@magmafunc($geqrf))
            magma_init()
            for i = 1:2                # 
                func(m, n, A, lda, tau, lwork, info)
                if i == 1
                    lwork = ceil(real(work[1]))
                    resize!(work, lwork)
                end
            end
            magma_finalize()
            A, tau
        end

    end

end
