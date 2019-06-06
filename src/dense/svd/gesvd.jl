include("magmaCgesvd.jl")
include("magmaDgesvd.jl")
include("magmaSgesvd.jl")
include("magmaZgesvd.jl")


# wrappers of the low-level MAGMA functions
for (function_name, element_type, singular_value_type) in
	((:magmaCgesvd,:ComplexF32,:Cfloat),
                      (:magmaZgesvd,:ComplexF32,:Cdouble))
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
		function gesvd!(jobu::Cint,
                      jobvt::Cint,
                      A::CuMatrix{$element_type},
                      ldu::Cint,
                      ldvt::Cint,
					  lwork::Cint,
					  rwork::PtrOrCuPtr{$singular_value_type})
			convert(Cint,jobu)
			convert(Cint,jobvt)
           m, n = size(A)
           # there should be some certain checking
		   # for the input matrix and arrays
           # before we can calculate further input
		   # and call the lower C function wrappers
           lda = max(1,stride(A,2))
           s = zeros(min(m,n))

		   # there should be some conditions for the size of U and VT
		   # but we will deal with them later
		   U = zeros($element_type,ldu,m)
		   VT = zeros($element_type,ldvt,n)

		   work = zeros($element_type,max(1,lwork))
		   info = zeros(Cint,1)

           $function_name(jobu, jobvt, m, n, A, lda, s, U, ldu, VT,
		    ldvt, work, lwork, rwork, info)
           return U, s, VT, work, info
	   end
    end
end

for (function_name, element_type, singular_value_type) in
	 ((:magmaDgesvd,:Cdouble,:Cdouble),
                      (:magmaSgesvd,:Cfloat,:Cfloat))
	@eval begin
		# magma_int_t magma_dgesvd 	( 	magma_vec_t  	jobu,
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
		# 	magma_int_t *  	info
		# )
		#Base.unsafe_convert(::CuPtr{$element_type},x::CuArray{$element_type})=Base.unsafe_convert(CuPtr{$element_type},Base.cconvert(CuPtr{($element_type)}x))
		#Base.unsafe_convert(::PtrOrCuPtr{$element_type},x::CuArray{$element_type})=Base.unsafe_convert(PtrOrCuPtr{$element_type},Base.cconvert(PtrOrCuPtr{($element_type)}x))
		function gesvd!(jobu,
                      jobvt,
                      A::Array{$element_type, 2},
                      ldu,
                      ldvt,
					  lwork)
			convert(Cint,jobu)
			convert(Cint,jobvt)
           m, n = size(A)
           # there should be some certain checking
		   # for the input matrix and arrays
           # before we can calculate further input
		   # and call the lower C function wrappers
           lda = max(1,stride(A,2))
           s = (zeros(min(m,n)))

		   # there should be some conditions for the size of U and VT
		   # but we will deal with them later
		   U = (zeros($element_type,ldu,m))
		   VT = (zeros($element_type,ldvt,n))

		   work = (zeros($element_type,max(1,lwork)))

		   info = (zeros(Cint,1))

		   $function_name(jobu, jobvt, m, n, A, lda, s, U, ldu, VT,
		    ldvt, work, lwork, info)
           return U, s, VT, work, info
	   end
    end
end
