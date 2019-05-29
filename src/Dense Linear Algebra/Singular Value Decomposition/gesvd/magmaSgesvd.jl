# magma_sgesvd wrapper
function magmaSgesvd(jobu, jobvt, m, n, A, lda, s, U,
	 ldu, VT, ldvt, work, lwork, info)
    ccall((:magma_sgesvd_gpu, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, PtrOrCuPtr{Float32},
			    Cint, PtrOrCuPtr{Float32}, PtrOrCuPtr{Float32},
				Cint, PtrOrCuPtr{Float32}, Cint, PtrOrCuPtr{Float32},
				Cint, PtrOrCuPtr{Cint}),# test: Ptr->PtrOrCuPtr

                jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt,
				work, lwork, info)
end
