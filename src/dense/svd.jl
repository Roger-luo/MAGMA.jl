for (gesvd, gesdd, elty, relty) in    ((:sgesvd, :sgesdd, :Float32, :Float32),
                            (:dgesvd, :dgesdd, :Float64, :Float64)),
                            interface in (:Matrix, :CuMatrix)
    @eval begin
    # magma_int_t magma_cgesvd 	( 	magma_vec_t  	jobu,
    # 	magma_vec_t  	jobvt,
    # 	magma_int_t  	m,
    # 	magma_int_t  	n,
    # 	magmaFloatComplex *  	A,
    # 	magma_int_t  	lda,
    # 	float *  	s,
    # 	magmaFloatComplex *  	U,
    # 	magma_int_t  	ldu,
    # 	magmaFloatComplex *  	VT,
    # 	magma_int_t  	ldvt,
    # 	magmaFloatComplex *  	work,
    # 	magma_int_t  	lwork,
    # 	float *  	rwork,
    # 	magma_int_t *  	info
    # )
        function magma_gesvd!(jobu::AbstractChar, jobvt::AbstractChar, A::$interface{$elty})

                # if A is not stored in the CPU, convert it to the CPU,
                # using Julia's native converter
                # instead of MAGMA's
                A = Matrix{$elty}(A)

                m, n    = size(A)
                minmn   = min(m, n)
                lda     = max(1, stride(A, 2))

                S       = similar(A, $relty, minmn)

                U       = similar(A, $elty, jobu == 'A' ? (m, m) : (jobu == 'S' ? (m, minmn) : (m, 0)))
                ldu     = max(1, stride(U, 2))

                VT      = similar(A, $elty, jobvt == 'A' ? (n, n) : (jobvt == 'S' ? (minmn, n) : (n, 0)))
                ldvt    = max(1, stride(VT, 2))

                work    = Vector{$elty}(undef, 1)

                lwork   = -1
                info    = Ref{Int}()

                jobu_magma      = char_to_magmaInt(jobu)
                jobvt_magma     = char_to_magmaInt(jobvt)

                for i in 1:2
                    # println("@magmafunc(gesvd): ", @magmafunc($gesvd))
                    func = eval(@magmafunc($gesvd))
                    func(jobu_magma, jobvt_magma,
                        m, n,
                        A, lda,
                        S,
                        U, ldu,
                        VT, ldvt,
                        work, lwork,
                        info)
                    # ccall((@magmafunc($gesvd), libmagma), Int,
                    #         (Int, Int,
                    #         Int, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{$relty},
                    #         Ptr{$elty}, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{Int}),
                    #         jobu_magma, jobvt_magma,
                    #         m, n,
                    #         A, lda,
                    #         S,
                    #         U, ldu,
                    #         VT, ldvt,
                    #         work, lwork,
                    #         info)
                    if i == 1
                        lwork = ceil(Int, real(work[1]))
                        resize!(work, lwork)
                    end
                end

                if jobu == 'O'
                    return (A, S, VT)
                elseif jobvt == 'O'
                    return (U, S, A)
                else
                    return (U, S, VT)
                end
        end

        # magma_int_t magma_cgesdd 	( 	magma_vec_t  	jobu,
        # 	magma_vec_t  	jobvt,
        # 	magma_int_t  	m,
        # 	magma_int_t  	n,
        # 	magmaFloatComplex *  	A,
        # 	magma_int_t  	lda,
        # 	float *  	s,
        # 	magmaFloatComplex *  	U,
        # 	magma_int_t  	ldu,
        # 	magmaFloatComplex *  	VT,
        # 	magma_int_t  	ldvt,
        # 	magmaFloatComplex *  	work,
        # 	magma_int_t  	lwork,
        # 	float *  	rwork,
        # 	magma_int_t *  	info
        # )
        function magma_gesdd!(job::AbstractChar, A::$interface{$elty})

                A = Matrix{$elty}(A)

                m, n    = size(A)
                minmn   = min(m, n)
                lda     = max(1, stride(A, 2))

                if job == 'A'
                    U  = similar(A, $elty, (m, m))
                    VT = similar(A, $elty, (n, n))
                elseif job == 'S'
                    U  = similar(A, $elty, (m, minmn))
                    VT = similar(A, $elty, (minmn, n))
                elseif job == 'O'
                    U  = similar(A, $elty, (m, m >= n ? 0 : m))
                    VT = similar(A, $elty, (n, m >= n ? n : 0))
                else
                    U  = similar(A, $elty, (m, 0))
                    VT = similar(A, $elty, (n, 0))
                end

                ldu = max(1, stride(U, 2))
                ldvt = max(1, stride(VT, 2))

                S       = similar(A, $relty, minmn)


                work    = Vector{$elty}(undef, 1)
                lwork   = Int(-1)
                iwork = Vector{Int}(undef, 8*minmn)
                info    = Ref{Int}()

                job_magma      = char_to_magmaInt(job)

                for i in 1:2
                    func = eval(@magmafunc($gesdd))
                    func(
                        job_magma,
                        m, n,
                        A, lda,
                        S,
                        U, ldu,
                        VT, ldvt,
                        work, lwork,
                        iwork, info
                    )
                    # ccall((@magmafunc($gesdd), libmagma), Int,
                    #         (Int,
                    #         Int, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{$relty},
                    #         Ptr{$elty}, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{$elty}, Int,
                    #         Ptr{Int}, Ptr{Int}),

                    #         job_magma,
                    #         m, n,
                    #         A, lda,
                    #         S,
                    #         U, ldu,
                    #         VT, ldvt,
                    #         work, lwork,
                    #         iwork, info)
                    if i == 1
                        lwork = ceil(Int, nextfloat(real(work[1])))
                        resize!(work, lwork)
                    end
                end

                if job == 'O'
                    if m >= n
                        return (A, S, VT)
                    else
                        # ()__
                        # ||::Z__
                        # ||::|:::Z____
                        # ||::|:::|====|
                        # ||==|===|====|
                        # ||""|===|====|
                        # ||  `"""|====|
                        # ||      `""""`
                        return (U, S, A)
                    end
                end
            return (U, S, VT)
        end
    end
end

for (gesvd, gesdd, elty, relty) in    ((:cgesvd, :cgesdd, :ComplexF32, :Float32),
                            (:zgesvd, :zgesdd, :ComplexF64, :Float64)),
                            interface in (:Matrix, :CuMatrix)
    @eval begin
    # magma_int_t magma_cgesvd 	( 	magma_vec_t  	jobu,
    # 	magma_vec_t  	jobvt,
    # 	magma_int_t  	m,
    # 	magma_int_t  	n,
    # 	magmaFloatComplex *  	A,
    # 	magma_int_t  	lda,
    # 	float *  	s,
    # 	magmaFloatComplex *  	U,
    # 	magma_int_t  	ldu,
    # 	magmaFloatComplex *  	VT,
    # 	magma_int_t  	ldvt,
    # 	magmaFloatComplex *  	work,
    # 	magma_int_t  	lwork,
    # 	float *  	rwork,
    # 	magma_int_t *  	info
    # )
        function magma_gesvd!(jobu::AbstractChar, jobvt::AbstractChar, A::$interface{$elty})

                # if A is not stored in the CPU, convert it to the CPU,
                # using Julia's native converter
                # instead of MAGMA's
                A = Matrix{$elty}(A)

                m, n    = size(A)
                minmn   = min(m, n)
                lda     = max(1, stride(A, 2))

                S       = similar(A, $relty, minmn)

                U       = similar(A, $elty, jobu == 'A' ? (m, m) : (jobu == 'S' ? (m, minmn) : (m, 0)))
                ldu     = max(1, stride(U, 2))

                VT      = similar(A, $elty, jobvt == 'A' ? (n, n) : (jobvt == 'S' ? (minmn, n) : (n, 0)))
                ldvt    = max(1, stride(VT, 2))

                work    = Vector{$elty}(undef, 1)

                rwork = Vector{$relty}(undef, 5minmn)
                lwork   = -1
                info    = Ref{Int}()

                jobu_magma      = char_to_magmaInt(jobu)
                jobvt_magma     = char_to_magmaInt(jobvt)

                for i in 1:2
                    func = eval(@magmafunc($gesvd))
                    func(
                        jobu_magma, jobvt_magma,
                        m, n,
                        A, lda,
                        S,
                        U, ldu,
                        VT, ldvt,
                        work, lwork,

                            rwork,

                        info
                    )
                    # ccall((@magmafunc($gesvd), libmagma), Int,
                    #     (Int, Int,
                    #     Int, Int,
                    #     Ptr{$elty}, Int,
                    #     Ptr{$relty},
                    #     Ptr{$elty}, Int,
                    #     Ptr{$elty}, Int,
                    #     Ptr{$elty}, Int, Ptr{$relty},
                    #     Ptr{Int}),
                    #     jobu_magma, jobvt_magma,
                    #     m, n,
                    #     A, lda,
                    #     S,
                    #     U, ldu,
                    #     VT, ldvt,
                    #     work, lwork,

                    #         rwork,

                    #     info)
                    if i == 1
                        lwork = ceil(Int, real(work[1]))
                        resize!(work, lwork)
                    end
                end

                if jobu == 'O'
                    return (A, S, VT)
                elseif jobvt == 'O'
                    return (U, S, A)
                else
                    return (U, S, VT)
                end
        end

        # magma_int_t magma_cgesdd 	( 	magma_vec_t  	jobu,
        # 	magma_vec_t  	jobvt,
        # 	magma_int_t  	m,
        # 	magma_int_t  	n,
        # 	magmaFloatComplex *  	A,
        # 	magma_int_t  	lda,
        # 	float *  	s,
        # 	magmaFloatComplex *  	U,
        # 	magma_int_t  	ldu,
        # 	magmaFloatComplex *  	VT,
        # 	magma_int_t  	ldvt,
        # 	magmaFloatComplex *  	work,
        # 	magma_int_t  	lwork,
        # 	float *  	rwork,
        # 	magma_int_t *  	info
        # )
        function magma_gesdd!(job::AbstractChar, A::$interface{$elty})

            A = Matrix{$elty}(A)

            m, n    = size(A)
            minmn   = min(m, n)
            lda     = max(1, stride(A, 2))

            if job == 'A'
                U  = similar(A, $elty, (m, m))
                VT = similar(A, $elty, (n, n))
            elseif job == 'S'
                U  = similar(A, $elty, (m, minmn))
                VT = similar(A, $elty, (minmn, n))
            elseif job == 'O'
                U  = similar(A, $elty, (m, m >= n ? 0 : m))
                VT = similar(A, $elty, (n, m >= n ? n : 0))
            else
                U  = similar(A, $elty, (m, 0))
                VT = similar(A, $elty, (n, 0))
            end

            ldu = max(1, stride(U, 2))
            ldvt = max(1, stride(VT, 2))

            S       = similar(A, $relty, minmn)


            work    = Vector{$elty}(undef, 1)
            lwork   = Int(-1)
            rwork = Vector{$relty}(undef, job == 'N' ?
                7*minmn : minmn*max(5*minmn+7, 2*max(m,n)+2*minmn+1))
            iwork = Vector{Int}(undef, 8*minmn)
            info    = Ref{Int}()

            job_magma      = char_to_magmaInt(job)

            for i in 1:2
                func = eval(@magmafunc($gesdd))
                func(
                    job_magma,
                    m, n,
                    A, lda,
                    S,
                    U, ldu,
                    VT, ldvt,
                    work, lwork,
                    rwork,
                    iwork, info
                )
                # ccall((@magmafunc($gesdd), libmagma), Int,
                #         (Int,
                #         Int, Int,
                #         Ptr{$elty}, Int,
                #         Ptr{$relty},
                #         Ptr{$elty}, Int,
                #         Ptr{$elty}, Int,
                #         Ptr{$elty}, Int,
                #         Ptr{$relty},
                #         Ptr{Int}, Ptr{Int}),

                #         job_magma,
                #         m, n,
                #         A, lda,
                #         S,
                #         U, ldu,
                #         VT, ldvt,
                #         work, lwork,
                #         rwork,
                #         iwork, info)
                if i == 1
                    lwork = ceil(Int, nextfloat(real(work[1])))
                    resize!(work, lwork)
                end
            end

            if job == 'O'
                if m >= n
                    return (A, S, VT)
                else
                    # ()__
                    # ||::Z__
                    # ||::|:::Z____
                    # ||::|:::|====|
                    # ||==|===|====|
                    # ||""|===|====|
                    # ||  `"""|====|
                    # ||      `""""`
                    return (U, S, A)
                end
            end
        return (U, S, VT)
        end
    end
end

for (gebrd, elty, relty) in    ((:sgebrd, :Float32, :Float32),
                            (:dgebrd, :Float64, :Float64),
                            (:cgebrd, :ComplexF32, :Float32),
                            (:zgebrd, :ComplexF64, :Float64)),
                            interface in (:Matrix, :CuMatrix)
    @eval begin
    # magma_int_t magma_cgebrd 	( 	magma_vec_t  	jobu,
    # 	magma_vec_t  	jobvt,
    # 	magma_int_t  	m,
    # 	magma_int_t  	n,
    # 	magmaFloatComplex *  	A,
    # 	magma_int_t  	lda,
    # 	float *  	s,
    # 	magmaFloatComplex *  	U,
    # 	magma_int_t  	ldu,
    # 	magmaFloatComplex *  	VT,
    # 	magma_int_t  	ldvt,
    # 	magmaFloatComplex *  	work,
    # 	magma_int_t  	lwork,
    # 	float *  	rwork,
    # 	magma_int_t *  	info
    # )
        function magma_gebrd!(A::$interface{$elty})

                A = Matrix{$elty}(A)

                m, n    = size(A)
                minmn   = min(m, n)
                lda     = max(1, stride(A, 2))

                d = similar(A, $relty, minmn)
                e = similar(A, $relty, minmn-1)

                tauq = similar(A, $elty, minmn)
                taup = similar(A, $elty, minmn)

                work    = Vector{$elty}(undef, 1)
                lwork   = Int(-1)
                info = Ref{Int}()

                for i in 1:2 # first call returns lwork as work[1]
                    func = eval(@magmafunc($gebrd))
                    func(
                        m, n, A, lda,
                        d, e, tauq, taup,
                        work, lwork, info
                    )
                    # ccall((@magmafunc($gebrd), libmagma), Int,
                    #     (Int, Int, Ptr{$elty}, Int,
                    #     Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ptr{$elty},
                    #     Ptr{$elty}, Int, Ptr{Int}),

                    #     m, n, A, lda,
                    #     d, e, tauq, taup,
                    #     work, lwork, info)
                    if i == 1 
                        lwork = ceil(Int, nextfloat(real(work[1])))
                        resize!(work, lwork)
                    end
                end

            return (A, d, e, tauq, taup)
        end
    end
end
