#import MAGMA
using MAGMA
using MAGMA.Utilities

#using MAGMA:gesvd!
using MAGMA: MagmaAllVec, gesvd!, libmagma



### some JuliaGPU packages, maybe useful (who knows)
using CUDAdrv
using CUDAapi
using CUDAnative
using CuArrays

#
using Test, LinearAlgebra

matrixToTest = rand(Float64, 2, 2)

right_answer = svd(matrixToTest).S
S = right_answer

jobu = MagmaAllVec
jobvt = MagmaAllVec

ldu=2
ldvt=2
lwork=134
success=magmaInit()

U, s, VT, work, info = gesvd!(jobu,jobvt,matrixToTest,ldu,ldvt,lwork)

diff = S .- s
error_value = norm(diff)

@test error_value < 1e-7
