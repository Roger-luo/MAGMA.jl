# magma_sgesvd wrapper
function magmaSgesvd(jobu, jobvt,
	 m, n,
	 A, lda,
	 s,
	 U, ldu,
	 VT, ldvt,
	 work, lwork, info)
    ccall((:magma_sgesvd, libmagma),
               Cint,
               (Cint, Cint, 		# jobu , jobvt
			    Cint, Cint, 		# m, n
				PtrOrCuPtr{Cfloat}, Cint, 	# A, lda
				PtrOrCuPtr{Cfloat}, 		# s
				PtrOrCuPtr{Cfloat}, Cint, 	# U, ldu
				PtrOrCuPtr{Cfloat}, Cint, 	# VT, ldvt
				PtrOrCuPtr{Cfloat}, Cint, PtrOrCuPtr{Cint}),# work, lwork, info

                jobu, jobvt,
				m, n,
				A, lda,
				s,
				U, ldu,
				VT, ldvt,
				work, lwork, info)
end
