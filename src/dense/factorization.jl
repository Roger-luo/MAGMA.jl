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

        function magma_geqrf!(A::Matrix{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m, n  = size(A)
            if length(tau) != min(m,n)
                throw(DimensionMismatch("tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2                # first call returns lwork as work[1]
                ccall((@blasfunc($geqrf), liblapack), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, A, max(1,stride(A,2)), tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

    end

end
