# magma_sgesvd wrapper
function magmaSgesvd(jobu, jobvt, m, n, A, lda, s, U,
	 ldu, VT, ldvt, work, lwork, info)
    ccall((:magma_sgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, PtrOrCuPtr{Cfloat},
			    Cint, PtrOrCuPtr{Cfloat}, PtrOrCuPtr{Cfloat},
				Cint, PtrOrCuPtr{Cfloat}, Cint, PtrOrCuPtr{Cfloat},
				Cint, PtrOrCuPtr{Cint}),# test: Ptr->PtrOrCuPtr

                jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt,
				work, lwork, info)
end
