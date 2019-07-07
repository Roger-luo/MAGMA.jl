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


for (fname, elty, relty) in    ((:sgesvd, :Float32, :Float32),
                            (:dgesvd, :Float64, :Float64),
                            (:cgesvd, :ComplexF32, :Float32),
                            (:zgesvd, :ComplexF64, :Float64))
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
        function gesvd!(jobu::AbstractChar, jobvt::AbstractChar, A::AbstractMatrix{$elty})
                m, n    = size(A)
                minmn   = min(m, n)
                lda     = max(1, stride(A, 2))

                S       = similar(A, $relty, minmn)

                U       = similar(A, $elty, jobu == 'A' ? (m, m) : (jobu == 'S' ? (m, minmn) : (m, 0)))
                ldu     = max(1, stride(U, 2))

                VT      = similar(A, $elty, jobvt == 'A' ? (n, n) : (jobvt == 'S' ? (minmn, n) : (n, 0)))
                ldvt    = max(1, stride(VT, 2))

                work    = Vector{$elty}(undef, 1)
                cmplx   = eltype(A) <: Complex
                if cmplx
                    rwork = Vector{$relty}(undef, 5minmn)
                end
                lwork   = -1
                info    = Ref{Cint}()

                jobu_magma      = char_to_magmaInt(jobu)
                jobvt_magma     = char_to_magmaInt(jobvt)

                for i in 1:2
                    if cmplx
                        ccall((@magmafunc($fname), libmagma), Cint,
                            (Cint, Cint,
                            Cint, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$relty},
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint, Ptr{$relty},
                            Ptr{Cint}),
                            jobu_magma, jobvt_magma,
                            m, n,
                            A, lda,
                            S,
                            U, ldu,
                            VT, ldvt,
                            work, lwork,

                                rwork,

                            info)
                    else
                        ccall((@magmafunc($fname), libmagma), Cint,
                            (Cint, Cint,
                            Cint, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$relty},
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{Cint}),
                            jobu_magma, jobvt_magma,
                            m, n,
                            A, lda,
                            S,
                            U, ldu,
                            VT, ldvt,
                            work, lwork,
                            info)
                    end
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
    end
end

for (fname, elty, relty) in    ((:sgesdd, :Float32, :Float32),
                            (:dgesdd, :Float64, :Float64),
                            (:cgesdd, :ComplexF32, :Float32),
                            (:zgesdd, :ComplexF64, :Float64))
    @eval begin
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
        function gesdd!(job::AbstractChar, A::AbstractMatrix{$elty})
                if isa(A, CuMatrix)
                    A = Matrix{$elty}(A)
                end

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
                cmplx   = eltype(A) <: Complex
                if cmplx
                    rwork = Vector{$relty}(undef, job == 'N' ? 7*minmn : minmn*max(5*minmn+7, 2*max(m,n)+2*minmn+1))
                end

                iwork = Vector{Cint}(undef, 8*minmn)
                info    = Ref{Cint}()

                job_magma      = char_to_magmaInt(job)

                for i in 1:2
                    if cmplx
                        ccall((@magmafunc($fname), libmagma), Cint,
                            (Cint,
                            Cint, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$relty},
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$relty},
                            Ptr{Cint}, Ptr{Cint}),

                            job_magma,
                            m, n,
                            A, lda,
                            S,
                            U, ldu,
                            VT, ldvt,
                            work, lwork,

                                rwork,

                            iwork, info)
                    else
                        ccall((@magmafunc($fname), libmagma), Cint,
                            (Cint,
                            Cint, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$relty},
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{$elty}, Cint,
                            Ptr{Cint}, Ptr{Cint}),

                            job_magma,
                            m, n,
                            A, lda,
                            S,
                            U, ldu,
                            VT, ldvt,
                            work, lwork,
                            iwork, info)
                    end
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
