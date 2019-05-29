# magma_cgesvd wrapper
function magmaCgesvd(jobu, jobvt, m, n, A, lda, s, U,
	ldu, VT, ldvt, work, lwork, rwork, info)
    ccall((:magma_cgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, PtrOrCuPtr{ComplexF32},
			    Cint, PtrOrCuPtr{Float32}, PtrOrCuPtr{ComplexF32},
				Cint, PtrOrCuPtr{ComplexF32}, Cint,
				PtrOrCuPtr{ComplexF32}, Cint, PtrOrCuPtr{Float32},
				PtrOrCuPtr{Cint}),
                jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt,
			    work, lwork, rwork, info)
end
