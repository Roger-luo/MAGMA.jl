using TensorOperations
using QuadGK
using LinearAlgebra
using Einsum

function itebd(G_list, l_list, U, chi)
    d = size(G_list[1])[1]

    for ibond in [0 1]
        # println(ibond)
        ia = mod(ibond, 2) + 1
        ib = mod(ibond + 1, 2) + 1

        chi1 = size(G_list[ia])[2]
        chi3 = size(G_list[ib])[3]

        # Construct theta
        Gia = G_list[ia]
        Gib = G_list[ib]
        println("GIA ", Gia)
        println("l list = ", l_list)
        Dia = diagm(0=>l_list[ia])
        Dib = diagm(0=>l_list[ib])
        println("Dia: ", Dia)
        theta = zeros((size(Dib, 1), size(Gia, 1), size(Gib, 1), size(Dib, 2)))
        @tensor theta[i, j, k, l] = Dib[i, b] * Gia[j, b, d] * Dia[d, e] * Gib[k, e, g] * Dib[g, l]

        # Apply U 
        @tensor theta[i, j, k, l] = theta[i, a, b, j] * reshape(U, (d, d, d, d))[a, b, k, l]

        # SVD
        theta = reshape(transpose(theta, (3, 1, 4, 2)), (d * chi1, d * chi3))
        X, Y, Z = svd(theta)
        Z = Z.transpose()

        chi2 = min([sum(Y > 10. ^ (-10)), chi_max])

        # Truncate
        l_list[ia] = Y[1:chi2] / sqrt(sum(Y[1:chi2] ^ 2))

        X = reshape(X[:, 0:chi2], (d, chi1, chi2))

        Dibinv = diag(l_list[ib] ^ (-1))

        temp = zeros(size(Dibinv, 1), size(X, 1), size(X, 3))
        @tensor temp[i, j, k]=Dibinv[i, a]*X[j, a, k]
        G_list[ia] = permutedims(temp, (2, 1, 3))

        Z = permutedims(reshape(Z[:, 1:chi2], (d, chi3, chi2)), (0, 2, 1))
        temp2 = zeros(size(Z, 1), size(Z, 2), size(Dibinv, 2))
        @tensor temp2[i, j, k]=Z[i, j, a]*Dibinv[a, k]
        G_list[ib] = temp2

        println("G list: ", G_list)
        
        # TODO: complete the codes

    end
end

function site_expectation_value(G_list, l_list, O)
    # Expectation value for a site operator
    E = []
    for isite in 1:2
        # TODO: complete the codes
    end
    
    return E
    
end

function bond_expectation_value(G_list, l_list, O)
    # Expectation value for a site operator
    E = []
    for ibond in 1:2
        
    end

    return E
end

# Define the simulation parameters
chi_max   = 10
δ       = 0.01
N       = 1000
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

# Initial state: |0000>
Ga = zeros((d, 1, 1))
Ga[1, 1, 1] = 1.
Gb = zeros((d, 1, 1))
Gb[1, 1, 1] = 1.
G_list = [Ga, Gb]

la = zeros(1)
la[1] = 1.
lb = zeros(1)
lb[1] = 1.
l_list = [la lb]

for step in 1:4
    itebd(G_list, l_list, U, chi_max)
end
println("m       = ", mean(site_expectation_value(G_list, l_list, sz)))
println("E_itebd = ", mean(bond_expectation_value(G_list, l_list, H)))

# Get the exact energy
function f(k, g)
    -2 * sqrt(1 + g * g - 2 * g * cos(k)) / pi / 2.
end
function ff(k)
    return f(k, g)
end
E0_exact = quadgk(ff, 0, pi)
println("E_exact = ", E0_exact)
