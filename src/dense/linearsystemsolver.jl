using Base: has_offset_axes

## (GE) general matrices, solvers with factorization, solver and inverse
# for (gels, gesv, getrs, getri, elty) in
#     ((:dgels,:dgesv,:dgetrs,:dgetri,:Float64),
#      (:sgels,:sgesv,:sgetrs,:sgetri,:Float32),
#      (:zgels,:zgesv,:zgetrs,:zgetri,:ComplexF64),
#      (:cgels,:cgesv,:cgetrs,:cgetri,:ComplexF32))
const function_list = ("gels", "gesv", "getrs", "getri", "getrf")
for type in magmaTypeList
    # create the symbols for element types
    elty = Symbol(type)
    # generate the symbol variables for our wrappers
    for func_name in function_list
        @eval $(Symbol(func_name)) = (Symbol(magmaTypeDict[$type], $func_name))
    end
    
    @eval begin
        #      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,INFO)
        # *     .. Scalar Arguments ..
        #       CHARACTER          TRANS
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
        function magma_gels!(trans::magma_trans_t, A::CuArray{$elty}, B::CuArray{$elty})
            m, n  = size(A)
            info  = Ref{Cint}()
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)

            for i = 1:2  # first call returns lwork as work[1]
                func = eval(@magmafunc_gpu($gels))
                func(
                    MagmaNoTrans, m, n, size(B,2),
                    A, max(1,stride(A,2)), B, max(1,stride(B,2)),
                    work, lwork, info
                )
                # ccall((@magmafunc_gpu($gels), libmagma), Cint,
                #       (Cint, Cint, Cint, Cint,
                #        PtrOrCuPtr{$elty}, Cint, PtrOrCuPtr{$elty}, Cint,
                #        Ptr{$elty}, Cint, Ptr{Cint}),
                #       MagmaNoTrans, m, n, size(B,2),
                #       A, max(1,stride(A,2)), B, max(1,stride(B,2)),
                #       work, lwork, info)

                if i == 1
                    lwork = ceil(Int, real(work[1]))
                    resize!(work, lwork)
                end
            end

            k   = min(m, n)
            F   = m < n ? tril(A[1:k, 1:k]) : triu(A[1:k, 1:k])
            ssr = Vector{$elty}(undef, size(B, 2))
            for i = 1:size(B,2)
                x = zero($elty)
                for j = k+1:size(B,1)
                    x += abs2(B[j,i])
                end
                ssr[i] = x
            end
            F, subsetrows(B, B, k), ssr
        end
        function magma_gels!(trans::magma_trans_t, A::Array{$elty}, B::Array{$elty})
            m, n  = size(A)
            info  = Ref{Cint}()
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)

            for i = 1:2  # first call returns lwork as work[1]
                func = eval(@magmafunc($gels))
                func(
                    MagmaNoTrans, m, n, size(B,2),
                    A, max(1,stride(A,2)), B, max(1,stride(B,2)),
                    work, lwork, info
                )
                # ccall((@magmafunc($gels), libmagma), Cint,
                #       (Cint, Cint, Cint, Cint,
                #        PtrOrCuPtr{$elty}, Cint, PtrOrCuPtr{$elty}, Cint,
                #        Ptr{$elty}, Cint, Ptr{Cint}),
                #       MagmaNoTrans, m, n, size(B,2),
                #       A, max(1,stride(A,2)), B, max(1,stride(B,2)),
                #       work, lwork, info)

                if i == 1
                    lwork = ceil(Int, real(work[1]))
                    resize!(work, lwork)
                end
            end

            k   = min(m, n)
            F   = m < n ? tril(A[1:k, 1:k]) : triu(A[1:k, 1:k])
            ssr = Vector{$elty}(undef, size(B, 2))
            for i = 1:size(B,2)
                x = zero($elty)
                for j = k+1:size(B,1)
                    x += abs2(B[j,i])
                end
                ssr[i] = x
            end
            F, subsetrows(B, B, k), ssr
        end

        #     Parameters
        # [in]	     n	    INTEGER The order of the matrix A. N >= 0.
        # [in]	     nrhs	INTEGER The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
        # [in,out]	 A	    COMPLEX_16 array, dimension (LDA,N). On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored.
        # [in]	     lda	INTEGER The leading dimension of the array A. LDA >= max(1,N).
        # [out]	     ipiv	INTEGER array, dimension (min(M,N)) The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).
        # [in,out]	 B	    COMPLEX_16 array, dimension (LDB,NRHS) On entry, the right hand side matrix B. On exit, the solution matrix X.
        # [in]	     ldb	INTEGER The leading dimension of the array B. LDB >= max(1,N).
        # [out]	     info	INTEGER
        #
        #               = 0: successful exit
        #               < 0: if INFO = -i, the i-th argument had an illegal value
        function magma_gesv!(A::CuArray{$elty}, B::CuArray{$elty})
            # require_one_based_indexing(A, B)
            # chkstride1(A, B)
            n = checksquare(A)
            isGPU = (A isa CuArray && B isa CuArray)
            # if size(B,1) != n
            #     throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
            # end
            lda  = max(1,stride(A,2))
            ldb  = max(1,stride(B,2))
            ipiv = similar(Matrix(A), Int32, n)
            info = Ref{Cint}()

            func = eval(@magmafunc_gpu($gesv))
            func(
                n, size(B,2),
                A, lda, ipiv,
                B, ldb, info
            )
            # ccall((@magmafunc_gpu($gesv), libmagma), Cint,
            #       (Cint, Cint,
            #        PtrOrCuPtr{$elty}, Cint, Ptr{Cint},
            #        PtrOrCuPtr{$elty}, Cint, Ptr{Cint}),
            #        n, size(B,2),
            #        A, lda, ipiv,
            #        B, ldb, info)
            # chkmagmaerror(info[])
            B, A, ipiv
        end
    
        function magma_gesv!(A::Array{$elty}, B::Array{$elty})
            # require_one_based_indexing(A, B)
            # chkstride1(A, B)
            n = checksquare(A)
            isGPU = (A isa CuArray && B isa CuArray)
            # if size(B,1) != n
            #     throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
            # end
            lda  = max(1,stride(A,2))
            ldb  = max(1,stride(B,2))
            ipiv = similar(Matrix(A), Int32, n)
            info = Ref{Cint}()

            func = eval(@magmafunc($gesv))
            func(
                n, size(B,2),
                A, lda, ipiv,
                B, ldb, info
            )
            # ccall((@magmafunc($gesv), libmagma), Cint,
            #       (Cint, Cint,
            #        PtrOrCuPtr{$elty}, Cint, Ptr{Cint},
            #        PtrOrCuPtr{$elty}, Cint, Ptr{Cint}),
            #        n, size(B,2),
            #        A, lda, ipiv,
            #        B, ldb, info)

            # chkmagmaerror(info[])
            B, A, ipiv
        end

        #     SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          TRANS
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function magma_getrs!(trans::magma_trans_t, A::CuMatrix{$elty}, ipiv::Array{Cint}, B::CuArray{$elty})
            # require_one_based_indexing(A, ipiv, B)
            # chktrans(trans)
            # chkstride1(A, B, ipiv)
            @assert !has_offset_axes(A, ipiv, B)
            n = checksquare(A)
            if n != size(B, 1)
                throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
            end
            nrhs = size(B, 2)
            lda  = max(1,stride(A,2))
            ldb  = max(1,stride(B,2))
            info = Ref{Cint}()
            
            func = eval(@magmafunc_gpu($getrs))
            func(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            B
            # println("\nThe info variable is: ", info, "\n")
        end
        function magma_getrs!(trans::magma_trans_t, A::Matrix{$elty}, ipiv::AbstractVector{Cint}, B::Array{$elty})
            magma_getrs!(trans, cu(A), ipiv, cu(B))
        end

        #*# Parameters of getri
        #*# [in]	n	INTEGER The order of the matrix A. N >= 0.
        #*# [in,out]	dA	COMPLEX array on the GPU, dimension (LDDA,N) On entry, the factors L and U from the factorization A = P*L*U as computed by CGETRF_GPU. On exit, if INFO = 0, the inverse of the original matrix A.
        #*# [in]	ldda	INTEGER The leading dimension of the array A. LDDA >= max(1,N).
        #*# [in]	ipiv	INTEGER array, dimension (N) The pivot indices from CGETRF; for 1 <= i <= N, row i of the matrix was interchanged with row IPIV(i).
        #*# [out]	dwork	(workspace) COMPLEX array on the GPU, dimension (MAX(1,LWORK))
        #*# [in]	lwork	INTEGER The dimension of the array DWORK. LWORK >= N*NB, where NB is the optimal blocksize returned by magma_get_cgetri_nb(n). 
        #*# Unlike LAPACK, this version does not currently support a workspace query, because the workspace is on the GPU.
        #*# [out]	info	INTEGER
        #*#         = 0: successful exit
        #*#         < 0: if INFO = -i, the i-th argument had an illegal value
        #*#         > 0: if INFO = i, U(i,i) is exactly zero; the matrix is singular and its cannot be computed.
        function magma_getri!(A::CuMatrix{$elty}, ipiv::Array{Int})
            # require_one_based_indexing(A, ipiv)
            # chkstride1(A, ipiv)
            n = checksquare(A)
            if n != length(ipiv)
                throw(DimensionMismatch("ipiv has length $(length(ipiv)), but needs $n"))
            end
            lda = max(1,stride(A, 2))

            func_nb = eval(@magmafunc_nb($getri))
            nb    = func_nb(n)
            lwork = ceil(Int, real(n * nb))
            work  = cu(Vector{$elty}(undef, max(1, lwork)))
            info  = Ref{Cint}()
            # for i = 1:2  # first call returns lwork as work[1]
                # ccall((@magmafunc_gpu($getri), libmagma), Cvoid, # ! MAGMA has no native cpu interface for getri
                #       (Ref{Cint}, PtrOrCuPtr{$elty}, Ref{Cint}, PtrOrCuPtr{Cint},
                #        PtrOrCuPtr{$elty}, Ref{Cint}, PtrOrCuPtr{Cint}),
                #       n, A, lda, ipiv, work, lwork, info)
                # # chkmagmaerror(info[])
            func = eval(@magmafunc_gpu($getri))
            func(n, A, lda, ipiv, work, lwork, info)
            #     if i == 1
            #         lwork = ceil(Int, real(work[1]))
            #         resize!(work, lwork)
            #     end
            # end
            A
        end
        function magma_getri!(A::Matrix{$elty}, ipiv::Array{Int})
            magma_getri!(cu(A), ipiv)
        end

        function magma_getrf!(A::AbstractMatrix{$elty})
            @assert !has_offset_axes(A)
            chkstride1(A)
            m, n = size(A)
            lda = max(1, stride(A, 2))
            ipiv = similar(A, Cint, min(m, n))
            info = Ref{Cint}()
            func = eval(@magmafunc($getrf))
            func(m, n, A, lda, ipiv, info)
            A, ipiv
        end
    end
end
