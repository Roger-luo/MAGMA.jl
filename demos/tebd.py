import numpy as np
from scipy.linalg import expm
from scipy import integrate


def itebd(G_list, l_list, U, chi_max):
    " Updates the G and s matrices using U and the TEBD protocol "
    d = G_list[0].shape[0]

    for ibond in [0, 1]:
        ia = np.mod(ibond, 2)
        ib = np.mod(ibond + 1, 2)

        chi1 = G_list[ia].shape[1]
        chi3 = G_list[ib].shape[2]

        # Construct theta
        theta = np.tensordot(np.diag(l_list[ib]), G_list[ia], axes=(1, 1))
        theta = np.tensordot(theta, np.diag(l_list[ia], 0), axes=(2, 0))
        theta = np.tensordot(theta, G_list[ib], axes=(2, 1))
        theta = np.tensordot(theta, np.diag(l_list[ib], 0), axes=(3, 0))

        # Apply U
        theta = np.tensordot(theta, np.reshape(U, (d, d, d, d)), axes=([1, 2], [0, 1]))

        # SVD
        theta = np.reshape(np.transpose(theta, (2, 0, 3, 1)), (d * chi1, d * chi3))  # ip a jp b
        X, Y, Z = np.linalg.svd(theta)
        Z = Z.T
        chi2 = np.min([np.sum(Y > 10. ** (-10)), chi_max])

        # Truncate
        l_list[ia] = Y[0:chi2] / np.sqrt(sum(Y[0:chi2] ** 2))

        X = np.reshape(X[:, 0:chi2], (d, chi1, chi2))
        G_list[ia] = np.transpose(np.tensordot(np.diag(l_list[ib] ** (-1)), X, axes=(1, 1)), (1, 0, 2))

        Z = np.transpose(np.reshape(Z[:, 0:chi2], (d, chi3, chi2)), (0, 2, 1))
        G_list[ib] = np.tensordot(Z, np.diag(l_list[ib] ** (-1)), axes=(2, 0))


def site_expectation_value(G_list, l_list, O):
    " Expectation value for a site operator "
    E = []
    for isite in range(0, 2):
        theta = np.tensordot(np.diag(l_list[np.mod(isite - 1, 2)]), G_list[isite], axes=(1, 1))
        theta = np.tensordot(theta, np.diag(l_list[isite]), axes=(2, 0))
        theta_O = np.tensordot(theta, O, axes=(1, 0)).conj()
        E.append(np.squeeze(np.tensordot(theta_O, theta, axes=([0, 1, 2], [0, 2, 1]))).item())
    return (E)


def bond_expectation_value(G_list, l_list, O):
    " Expectation value for a site operator "
    E = []
    for ibond in range(0, 2):
        ia = np.mod(ibond, 2)
        ib = np.mod(ibond + 1, 2)
        theta = np.tensordot(np.diag(l_list[ib]), G_list[ia], axes=(1, 1))
        theta = np.tensordot(theta, np.diag(l_list[ia], 0), axes=(2, 0))
        theta = np.tensordot(theta, G_list[ib], axes=(2, 1))
        theta = np.tensordot(theta, np.diag(l_list[ib], 0), axes=(3, 0))
        theta_O = np.tensordot(theta, np.reshape(O, (d, d, d, d)), axes=([1, 2], [0, 1])).conj()
        E.append(np.squeeze(np.tensordot(theta_O, theta, axes=([0, 1, 2, 3], [0, 3, 1, 2]))).item())

    return (E)


######## Define the simulation parameter ######################
chi_max = 10
delta = 0.01
N = 1000;
d = 2;
g = 0.5

########### Define Ising Hamiltonian and get U ################
sx = np.array([[0., 1.], [1., 0.]])
sz = np.array([[1., 0.], [0., -1.]])

H = -np.kron(sz, sz) + g * np.kron(sx, np.eye(2, 2))
U = expm(-delta * H)

############### Initial state : |0000> ########################
Ga = np.zeros((d, 1, 1), dtype=float)
Ga[0, 0, 0] = 1.
Gb = np.zeros((d, 1, 1), dtype=float)
Gb[0, 0, 0] = 1.
G_list = [Ga, Gb]

la = np.zeros(1)
la[0] = 1.
lb = np.zeros(1)
lb[0] = 1.
l_list = [la, lb]

for step in range(1, N):
    itebd(G_list, l_list, U, chi_max)
print("m       = ", np.mean(site_expectation_value(G_list, l_list, sz)))
print("E_itebd = ", np.mean(bond_expectation_value(G_list, l_list, H)))

############### Get the exact energy #########################
f = lambda k, g: -2 * np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)) / np.pi / 2.
E0_exact = integrate.quad(f, 0, np.pi, args=(g,))[0]
print("E_exact = ", E0_exact)
