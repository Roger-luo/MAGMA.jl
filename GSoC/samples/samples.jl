# sample code for magma binding


# MAGMA enum constants
# the whole file will be stored as enums.jl just like in JuliaGPU/MAGMA.jl

# MAGMA constants indicating the vectors status as input/output for some functions
# For example, the gesvd functions will use MagmaNoVec, MagmaSomeVec, MagmaAllVec and MagmaOverwriteVec to indicate the
# strategies that will be applied to the SVD U matrix and VT matrix in A = U Σ V**T (for MagmaOverwriteVec it is going to overwrite A)
const MagmaNoVec         = 301
const MagmaVec           = 302
const MagmaIVec          = 303
const MagmaAllVec        = 304
const MagmaSomeVec       = 305
const MagmaOverwriteVec  = 306
const MagmaBacktransVec = 307



# low-level wrappers of the MAGMA library

# magma_cgesvd wrapper
function magmaCgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
  @check ccall((:magma_cgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, CuPtr{ComplexF32}, Cint, PtrOrCuPtr{Cfloat}, CuPtr{ComplexF32}, Cint, CuPtr{ComplexF32}, Cint, CuPtr{ComplexF32}, Cint, PtrOrCuPtr{Cfloat}, CuPtr{Cint}),
               jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
end

# magma_dgesvd wrapper
function magmaDgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
  @check ccall((:magma_dgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, CuPtr{Cdouble}, Cint, PtrOrCuPtr{Cdouble}, CuPtr{Cdouble}, Cint, CuPtr{Cdouble}, Cint, CuPtr{Cdouble}, Cint, CuPtr{Cint}),
               jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
end

# magma_sgesvd wrapper
function magmaSgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
  @check ccall((:magma_sgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, CuPtr{Cfloat}, Cint, PtrOrCuPtr{Cfloat}, CuPtr{ComplexF32}, Cint, CuPtr{ComplexF32}, Cint, CuPtr{ComplexF32}, Cint, CuPtr{Cint}),
               jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
end

# magma_cgesvd wrapper
function magmaZgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
  @check ccall((:magma_zgesvd, libmagma),
               Cint,
               (Cint, Cint, Cint, Cint, CuPtr{ComplexF64}, Cint, PtrOrCuPtr{Cdouble}, CuPtr{ComplexF64}, Cint, CuPtr{ComplexF64}, Cint, CuPtr{ComplexF64}, Cint, PtrOrCuPtr{Cdouble}, CuPtr{Cint}),
               jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
end


# wrappers of the low-level MAGMA functions
for (function_name, element_type, singular_value_type) in ((:magmaCgesvd,:ComplexF32,:Cfloat),
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
           cutransa = cublasop(transa)
           cutransb = cublasop(transb)
           m, n = size(A)
           # there should be some certain checking for the input matrix and arrays
           # before we can calculate further input and call the lower C function wrappers
           lda = max(1,stride(A,2))
           s = CuArray{$singular_value_type}(0,min(m,n))

		   # there should be some conditions for the size of U and VT
		   # but we will deal with them later
		   U = zeros($element_type,ldu,m)
		   VT = zeros($element_type,ldvt,n)

		   work = zeros($element_type,max(1,lwork))
		   info = zeros(Cint,1)

           $function_name(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
           return U, s, VT, work, info
	   end
    end
end
for (function_name, element_type, singular_value_type) in ((:magmaDgesvd,:Cdouble,:Cdouble),
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
		function gesvd!(jobu::Cint,
                      jobvt::Cint,
                      A::CuMatrix{$element_type},
                      ldu::Cint,
                      ldvt::Cint,
					  lwork::Cint)
           cutransa = cublasop(transa)
           cutransb = cublasop(transb)
           m, n = size(A)
           # there should be some certain checking for the input matrix and arrays
           # before we can calculate further input and call the lower C function wrappers
           lda = max(1,stride(A,2))
           s = CuArray{$singular_value_type}(0,min(m,n))

		   # there should be some conditions for the size of U and VT
		   # but we will deal with them later
		   U = zeros($element_type,ldu,m)
		   VT = zeros($element_type,ldvt,n)

		   work = zeros($element_type,max(1,lwork))
		   info = zeros(Cint,1)

           $function_name(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
           return U, s, VT, work, info
	   end
    end
end
