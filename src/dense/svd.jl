for (gesvd, gesdd, elty, relty) in    ((:sgesvd, :sgesdd, :Float32, :Float32),
                            (:dgesvd, :dgesdd, :Float64, :Float64)),
                            interface in (:Matrix, :CuMatrix)
    @eval begin
        # *CGESVD computes the singular value decomposition (SVD) of a complex M-by-N matrix A, optionally computing the left and/or right singular vectors.

        # The SVD is written
        
        # A = U * SIGMA * conjugate-transpose(V)
        # where SIGMA is an M-by-N matrix which is zero except for its min(m,n) diagonal elements, U is an M-by-M unitary matrix, and V is an N-by-N unitary matrix. The diagonal elements of SIGMA are the singular values of A; they are real and non-negative, and are returned in descending order. The first min(m,n) columns of U and V are the left and right singular vectors of A.
        
        # Note that the routine returns VT = V**H, not V.
        
        # Parameters
        # [in]	jobu	magma_vec_t Specifies options for computing all or part of the matrix U:
        # = MagmaAllVec: all M columns of U are returned in array U:
        # = MagmaSomeVec: the first min(m,n) columns of U (the left singular vectors) are returned in the array U;
        # = MagmaOverwriteVec: the first min(m,n) columns of U (the left singular vectors) are overwritten on the array A;
        # = MagmaNoVec: no columns of U (no left singular vectors) are computed.
        # [in]	jobvt	magma_vec_t Specifies options for computing all or part of the matrix V**H:
        # = MagmaAllVec: all N rows of V**H are returned in the array VT;
        # = MagmaSomeVec: the first min(m,n) rows of V**H (the right singular vectors) are returned in the array VT;
        # = MagmaOverwriteVec: the first min(m,n) rows of V**H (the right singular vectors) are overwritten on the array A;
        # = MagmaNoVec: no rows of V**H (no right singular vectors) are computed. 
        # JOBVT and JOBU cannot both be MagmaOverwriteVec.
        # [in]	m	INTEGER The number of rows of the input matrix A. M >= 0.
        # [in]	n	INTEGER The number of columns of the input matrix A. N >= 0.
        # [in,out]	A	COMPLEX array, dimension (LDA,N) On entry, the M-by-N matrix A. On exit,
        # if JOBU = MagmaOverwriteVec, A is overwritten with the first min(m,n) columns of U (the left singular vectors, stored columnwise);
        # if JOBVT = MagmaOverwriteVec, A is overwritten with the first min(m,n) rows of V**H (the right singular vectors, stored rowwise);
        # if JOBU != MagmaOverwriteVec and JOBVT != MagmaOverwriteVec, the contents of A are destroyed.
        # [in]	lda	INTEGER The leading dimension of the array A. LDA >= max(1,M).
        # [out]	s	REAL array, dimension (min(M,N)) The singular values of A, sorted so that S(i) >= S(i+1).
        # [out]	U	COMPLEX array, dimension (LDU,UCOL) (LDU,M) if JOBU = MagmaAllVec or (LDU,min(M,N)) if JOBU = MagmaSomeVec.
        # If JOBU = MagmaAllVec, U contains the M-by-M unitary matrix U;
        # if JOBU = MagmaSomeVec, U contains the first min(m,n) columns of U (the left singular vectors, stored columnwise);
        # if JOBU = MagmaNoVec or MagmaOverwriteVec, U is not referenced.
        # [in]	ldu	INTEGER The leading dimension of the array U. LDU >= 1; if JOBU = MagmaSomeVec or MagmaAllVec, LDU >= M.
        # [out]	VT	COMPLEX array, dimension (LDVT,N)
        # If JOBVT = MagmaAllVec, VT contains the N-by-N unitary matrix V**H;
        # if JOBVT = MagmaSomeVec, VT contains the first min(m,n) rows of V**H (the right singular vectors, stored rowwise);
        # if JOBVT = MagmaNoVec or MagmaOverwriteVec, VT is not referenced.
        # [in]	ldvt	INTEGER The leading dimension of the array VT. LDVT >= 1;
        # if JOBVT = MagmaAllVec, LDVT >= N;
        # if JOBVT = MagmaSomeVec, LDVT >= min(M,N).
        # [out]	work	(workspace) COMPLEX array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK[0] returns the required LWORK.
        # [in]	lwork	INTEGER The dimension of the array WORK. If lwork = -1, a workspace query is assumed. The optimal size for the WORK array is calculated and stored in WORK[0], and no other work except argument checking is performed. 
        # Let mx = max(M,N) and mn = min(M,N). The threshold for mx >> mn is currently mx >= 1.6*mn. For job: N=None, O=Overwrite, S=Some, A=All. Paths below assume M >= N; for N > M swap jobu and jobvt. 
        # Because of varying nb for different subroutines, formulas below are an upper bound. Querying gives an exact number. The optimal block size nb can be obtained through magma_get_dgesvd_nb(M,N). For many cases, there is a fast algorithm, and a slow algorithm that uses less workspace. Here are sizes for both cases. 
        # Optimal lwork (fast algorithm) for mx >> mn: Path 1: jobu=N, jobvt=any 2*mn + 2*mn*nb Path 2: jobu=O, jobvt=N mn*mn + 2*mn + 2*mn*nb or mn*mn + max(2*mn + 2*mn*nb, mx*mn) Path 3: jobu=O, jobvt=A,S mn*mn + 2*mn + 2*mn*nb or mn*mn + max(2*mn + 2*mn*nb, mx*mn) Path 4: jobu=S, jobvt=N mn*mn + 2*mn + 2*mn*nb Path 5: jobu=S, jobvt=O 2*mn*mn + 2*mn + 2*mn*nb Path 6: jobu=S, jobvt=A,S mn*mn + 2*mn + 2*mn*nb Path 7: jobu=A, jobvt=N mn*mn + max(2*mn + 2*mn*nb, mn + mx*nb) Path 8: jobu=A, jobvt=O 2*mn*mn + max(2*mn + 2*mn*nb, mn + mx*nb) Path 9: jobu=A, jobvt=A,S mn*mn + max(2*mn + 2*mn*nb, mn + mx*nb) for mx >= mn, but not mx >> mn: Path 10: jobu=any, jobvt=any 2*mn + (mx + mn)*nb 
        # Optimal lwork (slow algorithm) for mx >> mn: Path 1: jobu=N, jobvt=any n/a Path 2: jobu=O, jobvt=N 2*mn + (mx + mn)*nb Path 3-9: 2*mn + max(2*mn*nb, mx*nb) for mx >= mn, but not mx >> mn: Path 10: jobu=any, jobvt=any n/a 
        # MAGMA requires the optimal sizes above, while LAPACK has the same optimal sizes but the minimum sizes below. 
        # LAPACK minimum lwork (fast algorithm) for mx >> mn: Path 1: jobu=N, jobvt=any 3*mn Path 2: jobu=O, jobvt=N mn*mn + 3*mn Path 3: jobu=O, jobvt=A,S mn*mn + 3*mn Path 4: jobu=S, jobvt=N mn*mn + 3*mn Path 5: jobu=S, jobvt=O 2*mn*mn + 3*mn Path 6: jobu=S, jobvt=A,S mn*mn + 3*mn Path 7: jobu=A, jobvt=N mn*mn + max(3*mn, mn + mx) Path 8: jobu=A, jobvt=O 2*mn*mn + max(3*mn, mn + mx) Path 9: jobu=A, jobvt=A,S mn*mn + max(3*mn, mn + mx) for mx >= mn, but not mx >> mn: Path 10: jobu=any, jobvt=any 2*mn + mx 
        # LAPACK minimum lwork (slow algorithm) for mx >> mn: Path 1: jobu=N, jobvt=any n/a Path 2-9: 2*mn + mx for mx >= mn, but not mx >> mn: Path 10: jobu=any, jobvt=any n/a
        # rwork	(workspace) REAL array, dimension (5*min(M,N)) On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the unconverged superdiagonal elements of an upper bidiagonal matrix B whose diagonal is in S (not necessarily sorted). B satisfies A = U * B * VT, so it has the same singular values as A, and singular vectors related by U and VT.
        # [out]	info	INTEGER
        # = 0: successful exit.
        # < 0: if INFO = -i, the i-th argument had an illegal value.
        # > 0: if CBDSQR did not converge, INFO specifies how many superdiagonals of an intermediate bidiagonal form B did not converge to zero. See the description of RWORK above for details.
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
                info    = Ref{Cint}()

                jobu_magma      = char_to_magmaInt(jobu)
                jobvt_magma     = char_to_magmaInt(jobvt)
                func = eval(@magmafunc($gesvd))

                magma_init()
                for i in 1:2
                    func(jobu_magma, jobvt_magma,
                        m, n,
                        A, lda,
                        S,
                        U, ldu,
                        VT, ldvt,
                        work, lwork,
                        info)
                    # ccall((@magmafunc($gesvd), libmagma), Cint,
                    #         (Cint, Cint,
                    #         Cint, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$relty},
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{Cint}),
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
                magma_finalize()

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
                lwork   = Cint(-1)
                iwork = Vector{Cint}(undef, 8*minmn)
                info    = Ref{Cint}()

                job_magma      = char_to_magmaInt(job)
                func = eval(@magmafunc($gesdd))

                magma_init()
                for i in 1:2
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
                    # ccall((@magmafunc($gesdd), libmagma), Cint,
                    #         (Cint,
                    #         Cint, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$relty},
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{$elty}, Cint,
                    #         Ptr{Cint}, Ptr{Cint}),

                    #         job_magma,
                    #         m, n,
                    #         A, lda,
                    #         S,
                    #         U, ldu,
                    #         VT, ldvt,
                    #         work, lwork,
                    #         iwork, info)
                    if i == 1
                        lwork = ceil(Cint, nextfloat(real(work[1])))
                        resize!(work, lwork)
                    end
                end
                magma_finalize()

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
                info    = Ref{Cint}()

                jobu_magma      = char_to_magmaInt(jobu)
                jobvt_magma     = char_to_magmaInt(jobvt)

                func = eval(@magmafunc($gesvd))
                for i in 1:2
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
                    # ccall((@magmafunc($gesvd), libmagma), Cint,
                    #     (Cint, Cint,
                    #     Cint, Cint,
                    #     Ptr{$elty}, Cint,
                    #     Ptr{$relty},
                    #     Ptr{$elty}, Cint,
                    #     Ptr{$elty}, Cint,
                    #     Ptr{$elty}, Cint, Ptr{$relty},
                    #     Ptr{Cint}),
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
            lwork   = Cint(-1)
            rwork = Vector{$relty}(undef, job == 'N' ?
                7*minmn : minmn*max(5*minmn+7, 2*max(m,n)+2*minmn+1))
            iwork = Vector{Cint}(undef, 8*minmn)
            info    = Ref{Cint}()

            job_magma      = char_to_magmaInt(job)
            func = eval(@magmafunc($gesdd))

            for i in 1:2
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
                # ccall((@magmafunc($gesdd), libmagma), Cint,
                #         (Cint,
                #         Cint, Cint,
                #         Ptr{$elty}, Cint,
                #         Ptr{$relty},
                #         Ptr{$elty}, Cint,
                #         Ptr{$elty}, Cint,
                #         Ptr{$elty}, Cint,
                #         Ptr{$relty},
                #         Ptr{Cint}, Ptr{Cint}),

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
                    lwork = ceil(Cint, nextfloat(real(work[1])))
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
                lwork   = Cint(-1)
                info = Ref{Cint}()
                func = eval(@magmafunc($gebrd))

                for i in 1:2 # first call returns lwork as work[1]
                    func(
                        m, n, A, lda,
                        d, e, tauq, taup,
                        work, lwork, info
                    )
                    # ccall((@magmafunc($gebrd), libmagma), Cint,
                    #     (Cint, Cint, Ptr{$elty}, Cint,
                    #     Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ptr{$elty},
                    #     Ptr{$elty}, Cint, Ptr{Cint}),

                    #     m, n, A, lda,
                    #     d, e, tauq, taup,
                    #     work, lwork, info)
                    if i == 1 
                        lwork = ceil(Cint, nextfloat(real(work[1])))
                        resize!(work, lwork)
                    end
                end

            return (A, d, e, tauq, taup)
        end
    end
end
function magma_gebrd!(A::AbstractMatrix{T}) where T
    println("Begin to initialize\n")
    magma_init()
    magma_gebrd!(A)
    magma_finalize()
end
