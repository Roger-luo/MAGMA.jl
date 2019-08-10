## (GE) general matrices, solvers with factorization, solver and inverse
for (gels, gesv, getrs, getri, elty) in
    ((:dgels,:dgesv,:dgetrs,:dgetri,:Float64),
     (:sgels,:sgesv,:sgetrs,:sgetri,:Float32),
     (:zgels,:zgesv,:zgetrs,:zgetri,:ComplexF64),
     (:cgels,:cgesv,:cgetrs,:cgetri,:ComplexF32))
    @eval begin
        #      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,INFO)
        # *     .. Scalar Arguments ..
        #       CHARACTER          TRANS
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
        function magma_gels!(trans::Int, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            # require_one_based_indexing(A, B)
            #chktrans(trans)
            # chkstride1(A, B)
            btrn  = trans == 'T' # however currently the MAGMA only handles 'N'
            if btrn
                println("Now the MAGMA does not handle 'T' actually~")
            end
            m, n  = size(A)
            if size(B,1) != (btrn ? n : m)
                throw(DimensionMismatch("matrix A has dimensions ($m,$n), transposed: $btrn, but leading dimension of B is $(size(B,1))"))s
            end
            info  = Ref{Cint}()
            work  = Vector{$elty}(undef, 1)
            lwork = Cint(-1)
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@magmafunc($gels), libmagma), Cint,
                      (Cint, Cint, Cint, Cint,
                       Ptr{$elty}, Cint, Ptr{$elty}, Cint,
                       Ptr{$elty}, Cint, Ptr{Cint}),
                      MagmaNoTrans, m, n, size(B,2), A, max(1,stride(A,2)),
                      B, max(1,stride(B,2)), work, lwork, info)
                # chkmagmaerror(info[])
                if i == 1
                    lwork = Cint(real(work[1]))
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

        # SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function magma_gesv!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n = checksquare(A)
            if size(B,1) != n
                throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
            end
            ipiv = similar(A, Cint, n)
            info = Ref{Cint}()
            ccall((@magmafunc($gesv), libmagma), Cvoid,
                  (Ref{Cint}, Ref{Cint}, Ptr{$elty}, Ref{Cint}, Ptr{Cint},
                   Ptr{$elty}, Ref{Cint}, Ptr{Cint}),
                  n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info)
            chkmagmaerror(info[])
            B, A, ipiv
        end

        #     SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          TRANS
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function magma_getrs!(trans::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{Cint}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chktrans(trans)
            chkstride1(A, B, ipiv)
            n = checksquare(A)
            if n != size(B, 1)
                throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
            end
            nrhs = size(B, 2)
            info = Ref{Cint}()
            ccall((@magmafunc($getrs), libmagma), Cvoid,
                  (Ref{UInt8}, Ref{Cint}, Ref{Cint}, Ptr{$elty}, Ref{Cint},
                   Ptr{Cint}, Ptr{$elty}, Ref{Cint}, Ptr{Cint}),
                  trans, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info)
            chkmagmaerror(info[])
            B
        end

        #     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        #*     .. Scalar Arguments ..
        #      INTEGER            INFO, LDA, LWORK, N
        #*     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function magma_getri!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{Cint})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            if n != length(ipiv)
                throw(DimensionMismatch("ipiv has length $(length(ipiv)), but needs $n"))
            end
            lda = max(1,stride(A, 2))
            lwork = Cint(-1)
            work  = Vector{$elty}(undef, 1)
            info  = Ref{Cint}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@magmafunc($getri), libmagma), Cvoid,
                      (Ref{Cint}, Ptr{$elty}, Ref{Cint}, Ptr{Cint},
                       Ptr{$elty}, Ref{Cint}, Ptr{Cint}),
                      n, A, lda, ipiv, work, lwork, info)
                chkmagmaerror(info[])
                if i == 1
                    lwork = Cint(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A
        end
    end
end
