# magma_cgesvd wrapper
function magmaZgesvd(jobu, jobvt, m, n, A, lda, s, U,
	 ldu, VT, ldvt, work, lwork, rwork, info)
    ccall((:magma_zgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, PtrOrCuPtr{ComplexF64},
			    Cint, PtrOrCuPtr{Float64}, PtrOrCuPtr{ComplexF64},
				Cint, PtrOrCuPtr{ComplexF64}, Cint, PtrOrCuPtr{ComplexF64},
				Cint, PtrOrCuPtr{Float64}, PtrOrCuPtr{Cint}),
                jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt,
			    work, lwork, rwork, info)
end
