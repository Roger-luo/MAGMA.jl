const function_list_eigenvalues = ("geev", "gehrd", 
                                   "heevx", "heevd", "heevdx", "heevr", "hetrd",
                                   "hegvx", "hegvd", "hegvdx", "hegvr")
macro char_to_vec(op)
    return op == 'V' ? MagmaVec : MagmaNoVec
end
macro type_to_rtype(type)
    return type == Float64 || type == ComplexF64 ? Float64 : Float32
end

for type in magmaTypeList_real
    # create the symbols for element types
    elty = Symbol(type)
    relty= Symbol(@type_to_rtype(type))
    # generate the symbol variables for our wrappers
    for func_name in function_list_eigenvalues
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end

    @eval begin

        function magma_geev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty})
            magma_init()
            chkstride1(A)
            n = checksquare(A)
            chkfinite(A) # balancing routines don't support NaNs and Infs
            lvecs = jobvl == 'V'
            rvecs = jobvr == 'V'
            VL    = similar(A, $elty, (n, lvecs ? n : 0))
            VR    = similar(A, $elty, (n, rvecs ? n : 0))
            WR    = similar(A, $elty, n)
            WI    = similar(A, $elty, n)
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)
            info  = Ref{Cint}()
            lda   = max(1, stride(A, 2))

            jobvl = @char_to_vec(jobvl)
            jobvr = @char_to_vec(jobvr)
            func  = eval(@magmafunc($geev))
            for i = 1:2  # first call returns lwork as work[1]
                func(jobvl, jobvr, n, A, lda, WR, WI, VL, n,
                        VR, n, work, lwork, info)
                if i == 1
                    lwork = ceil(Int, real(work[1]))
                    resize!(work, lwork)
                end
            end
            magma_finalize()
            (WR, WI, VL, VR)
        end

    end

end

for type in magmaTypeList_complex
    # create the symbols for element types
    elty = Symbol(type)
    relty= Symbol(@type_to_rtype(type))
    # generate the symbol variables for our wrappers
    for func_name in function_list_eigenvalues
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end

    @eval begin

        function magma_geev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty})
            magma_init()
            chkstride1(A)
            n = checksquare(A)
            chkfinite(A) # balancing routines don't support NaNs and Infs
            lvecs = jobvl == 'V'
            rvecs = jobvr == 'V'
            VL    = similar(A, $elty, (n, lvecs ? n : 0))
            VR    = similar(A, $elty, (n, rvecs ? n : 0))
            W     = similar(A, $elty, n)
            rwork = similar(A, $relty, 2n)
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)
            info  = Ref{Cint}()
            lda   = max(1, stride(A, 2))

            jobvl = @char_to_vec(jobvl)
            jobvr = @char_to_vec(jobvr)
            func  = eval(@magmafunc($geev))
            for i = 1:2  # first call returns lwork as work[1]
                func(jobvl, jobvr, n, A, lda, W, VL, n, VR, n,
                    work, lwork, rwork, info)
                if i == 1
                    lwork = ceil(Int, real(work[1]))
                    resize!(work, lwork)
                end
            end
            magma_finalize()
            (W, VL, VR)
        end

    end

end

for type in magmaTypeList
    # create the symbols for element types
    elty = Symbol(type)
    relty= Symbol(@type_to_rtype(type))
    # generate the symbol variables for our wrappers
    for func_name in function_list_eigenvalues
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end

    @eval begin

        function magma_gehrd!(ilo::Integer, ihi::Integer, A::AbstractMatrix{$elty})
            magma_init()
            chkstride1(A)
            n = checksquare(A)
            chkfinite(A) # balancing routines don't support NaNs and Infs
            tau = similar(A, $elty, max(0,n - 1))
            work = Vector{$elty}(undef, 1)
            lwork = Cint(-1)
            info = Ref{Cint}()
            func  = eval(@magmafunc($gehrd))    
            nb    = eval(@magmafunc_nb($gehrd))(n)
            T     = CuArray{$elty}(undef, nb*n)
            for i = 1:2  # first call returns lwork as work[1]
                func(n, ilo, ihi, A,
                    max(1, stride(A, 2)), tau, work, lwork, T, 
                    info)
                if i == 1
                    lwork = ceil(Int, real(work[1]))
                    resize!(work, lwork)
                end
            end
            magma_finalize()
            A, tau
        end

    end

end
