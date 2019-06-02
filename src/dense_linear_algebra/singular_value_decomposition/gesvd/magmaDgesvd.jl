# magma_dgesvd wrapper
function magmaDgesvd(jobu, jobvt, m, n, A, lda, s, U,
	 ldu, VT, ldvt, work, lwork, info)
    ccall((:magma_dgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, PtrOrCuPtr{Cdouble},
			    Cint, PtrOrCuPtr{Cdouble}, PtrOrCuPtr{Cdouble},
				Cint, PtrOrCuPtr{Cdouble}, Cint,
				PtrOrCuPtr{Cdouble}, Cint, PtrOrCuPtr{Cint}),

                jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt,
			    work, lwork, info)
end
