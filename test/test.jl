#print(1+1)
print("Test begins!\n")

jobu = MagmaAllVec
jobvt = MagmaAllVec
A = cu(zeros(2,2))
A[1,1]=2.0
A[1,2]=1.0
A[2,1]=1.0
A[2,2]=2.0
ldu=2
ldvt=2
lwork=134
success=magmaInit()

print("Magma Initialization success=", success)
print('\n')

U, s, VT, work, info = gesvd!(jobu,jobvt,A,ldu,ldvt,lwork)

print("U=")
print(U)
print('\n')

print("s=")
print(s)
print('\n')

print("V**T=")
print(VT)
print('\n')

print("work=")
print(work)
print('\n')

print("info=")
print(info)
print('\n')

print("Test ends.")
