using TensorOperations
using QuadGK
using LinearAlgebra
using Einsum

function itebd(G_list, l_list, U, χ)
    d = size(G_list[1])[1]

    for ibond in [0 1]
        # println(ibond)
        ia = mod(ibond, 2) + 1
        ib = mod(ibond + 1, 2) + 1

        χ1 = size(G_list[ia])[2]
        χ3 = size(G_list[ib])[3]

        # Construct theta
        Gia = G_list[ia]
        Gib = G_list[ib]
        Dia = diagm(0=>l_list[ib])

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
χ_max   = 10
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
    itebd(G_list, l_list, U, χ_max)
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
