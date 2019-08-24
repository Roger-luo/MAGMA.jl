using TensorOperations
using QuadGK
using LinearAlgebra
using LinearAlgebra.LAPACK
using Einsum
using Statistics

function itebd(G_list, l_list, U, chi)
    d = size(G_list[1])[1]

    loop = 0
    for ibond in 1:2
        loop +=1
        # println("Loop number: ", loop)
        # println("Before loop G list: ", G_list)
        # println(ibond)
        ia = ibond
        ib = mod(ibond, 2) + 1
        # println("ia = ", ia, "\nib = ", ib)

        chi1 = size(G_list[ia])[2]
        chi3 = size(G_list[ib])[3]
        # println(chi1, chi3)

        # Construct theta
        Gia = G_list[ia]
        Gib = G_list[ib]
        # println("GIA ", Gia)
        # println("l list = ", l_list)
        Dia = diagm(0=>l_list[ia])
        Dib = diagm(0=>l_list[ib])
        # println("Dia: ", Dia)
        # theta0 = zeros((size(Dib, 1), size(Gia, 1), size(Gib, 1), size(Dib, 2)))
        # theta  = zeros((size(Dib, 1), size(Dib, 2), size(Gia, 1), size(Gib, 1)))
        # println("GList ia",Gia)
        # println("GList ib",Gib)
        # println("DList ia",Dia)
        # println("DList ib",Dib)
        
        # @tensor theta0[i, j, k, l] := Dib[i, b] * Gia[j, b, d] * Dia[d, e] * Gib[k, e, g] * Dib[g, l]
        # println("theta0 = ", theta0)
        # println("theta = ", theta)
        # println("Size of theta0 = ", size(theta0))

        # Apply U 
        # tempU = reshape(U, (d,d,d,d))
        # println("Size of temp U = ", size(tempU))
        # println("Size of theta = ", size(theta))
        # @tensor theta[i, j, k, l] := theta0[i, a, b, j] * tempU[a, b, k, l]
        
        """Multiple steps"""
        # @tensor theta[i, j, k] := Dib[i, a] * Gia[j, a, k]
        # # println("theta[i, j, k] := Dib[i, a] * Gia[j, a, k]\n", theta)
        # @tensor theta[i, j, k] := theta[i, j, a] * Dia[a, k]
        # # println("theta[i, j, k] := theta[i, j, a] * Dia[a, k]\n", theta)
        # @tensor theta[i, j, k, l] := theta[i, j, a] * Gib[k, a, l]
        # # println("theta[i, j, k, l] := theta[i, j, a] * Gib[k, a, l]\n", theta)
        # @tensor theta[i, j, k, l] := theta[i, j, k, a] * Dib[a, l]
        # # println("theta[i, j, k, l] := theta[i, j, k, a] * Dib[a, l]\n", theta)
        # Ur = reshape(U, (d,d,d,d))

        # # theta = permutedims(theta, (4, 3, 2, 1))
        # Ur = permutedims(Ur, (4, 3, 2, 1))
        # # println("Ur = ", Ur)
        # @tensor theta[i, j, k, l] := theta[i, a1, a2, j] * Ur[a1, a2, k, l]

        """One step"""
        @tensor theta[i, j, k, l] := Dib[i, b] * Gia[a, b, h] * Dia[h, e] * Gib[f, e, g] * Dib[g, j] * reshape(U, (d,d,d,d))[a, f, k, l]

        # SVD
        theta = reshape(permutedims(theta, (3, 1, 4, 2)), (d * chi1, d * chi3))
        
        # println("Theta reshaped = ", theta, " of the type", typeof(theta))
        # println("The size of theta is ", size(theta))
        # F = svd(theta)
        # X = F.U
        # Y = F.S
        # Z = F.Vt
        (X, Y, Z) = LAPACK.gesdd!('A', Matrix(theta))
        # println("X=",X)
        # println("Y=",Y)
        # println("Z=",Z)
        Z = permutedims(Z)
        notzero = Y.> (10^-10)
        # println("notezero=",notzero)
        # println("Sum = ", sum(notzero))
        chi2 = minimum([sum(notzero), chi_max])

        # Truncate
        l_list[ia] = Y[1:chi2] / sqrt(sum(Y[1:chi2] .^ 2))

        X = reshape(X[:, 1:chi2], (d, chi1, chi2))

        Dibinv = diagm(0=>l_list[ib] .^ (-1))
        @tensor temp[i, j, k]:=Dibinv[i, a]*X[j, a, k]
        G_list[ia] = permutedims(temp, (2, 1, 3))

        Z = permutedims(reshape(Z[:, 1:chi2], (d, chi3, chi2)), (1, 3, 2))
        temp2 = zeros(size(Z, 1), size(Z, 2), size(Dibinv, 2))
        @tensor temp2[i, j, k]=Z[i, j, a]*Dibinv[a, k]
        G_list[ib] = temp2

        # println("After loop G list: ", G_list)
        # println(" ")
        

    end
end

function site_expectation_value(G_list, l_list, O)
    # Expectation value for a site operator
    E = []
    for isite in 1:2
        Gisite = G_list[isite]
        Disite = diagm(0=>l_list[mod(isite, 2) + 1])
        Disite2= diagm(0=>l_list[isite])
        # theta  = zeros(size(Disite, 1), size(Gisite, 1), size(Disite2, 2))
        @tensor theta[i, j, k] := Disite[i, a]*Gisite[j, a, b]*Disite2[b, k]

        # theta_O= zeros(size(theta, 1), size(theta, 3), size(O, 2))
        @tensor theta_O[i, j, k] := theta[i, a, j] * O[a, k]
        conj!(theta_O)

        @tensor E_cal[]:=theta_O[i, j, k] * theta[i, k, j]
        push!(E, E_cal)
    end

    # println("E = ", E)
    
    return E
    
end

function bond_expectation_value(G_list, l_list, O)
    # Expectation value for a site operator
    E = []
    for ibond in 1:2
        ia = ibond
        ib = mod(ibond, 2)+1

        Gia= G_list[ia]
        Gib= G_list[ib]
        Dia= diagm(0=>l_list[ia])
        Dib= diagm(0=>l_list[ib])
        
        @tensor theta[i, j, k, l]:=Dib[i, a]*Gia[j, a, b]*Dia[b, d]*Gib[k, d, f]*Dib[f, l]
        O2 = reshape(O, (d,d,d,d))
        @tensor theta_O[i,j,k,l]:=theta[i, a1, a2, j]*O2[a1, a2, k, l]
        conj!(theta_O)
        @tensor E_cal[]:=theta_O[a0, a1, a2, a3]*theta[a0, a2, a3, a1]
        push!(E, E_cal)
    end

    return E
end

t0 = time()
# Define the simulation parameters
chi_max   = 10
δ       = 0.01
N       = 200
d       = 2
g       = 0.5

# Define Ising Hamiltonian and get using
id      = diagm(0=>ones(2))
sx      = [ [0.     1.  ]
            [1.     0.  ]]
sz      = [ [1.     0.  ]
            [0.    -1.  ]]

H       = -1 * kron(sz, sz) + g * kron(sx, id)
U       = exp(- δ * H )
# println("U=",U)
# Initial state: |0000>
Ga = zeros(Float64, (d, 1, 1))
Ga[1, 1, 1] = 1.
Gb = zeros(Float64, (d, 1, 1))
Gb[1, 1, 1] = 1.
G_list = [Ga, Gb]

la = zeros(1)
la[1] = 1.
lb = zeros(1)
lb[1] = 1.
l_list = [la, lb]

for step in 2:N
    @eval begin
    itebd(G_list, l_list, U, chi_max)
    end
end

println("m       = ", mean(site_expectation_value(G_list, l_list, sz)))
println("E_itebd = ", mean(bond_expectation_value(G_list, l_list, H)))

# Get the exact energy

function ff(k)
    return -2 * sqrt(1 + g * g - 2 * g * cos(k)) / pi / 2.
end
E0_exact = quadgk(ff, 0, pi)
println("E_exact = ", E0_exact)
# println("PROGRAMME FINISHED!!!\n\n\n")
t1 = time()
println("Time elapsed = ", t1-t0)
