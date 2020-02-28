import numpy as np
import itertools

from qutip import basis, tensor, Options, mesolve
from qutip import qeye, sigmax, sigmay, sigmaz


def Is(i):
    return [qeye(2) for j in range(0, i)]


def Sx(N, i):
    return tensor(Is(i) + [sigmax()] + Is(N - i - 1))


def Sy(N, i):
    return tensor(Is(i) + [sigmay()] + Is(N - i - 1))


def Sz(N, i):
    return tensor(Is(i) + [sigmaz()] + Is(N - i - 1))


def bloch(zero, one):
    z0 = zero
    z1 = one
    xp = (zero + one).unit()
    xm = (zero - one).unit()
    yp = (zero + 1j*one).unit()
    ym = (zero - 1j*one).unit()
    return z0, z1, xp, xm, yp, ym


def deriveUnitary(H, t):
    return (-1j*(t/2.)*H).expm()


def pmeasurement(psi, register, normalize=True):
    # get dimensions
    dim2 = psi.dims[1][0]
    # how many qubits get measured
    N = len(register)
    M = np.count_nonzero(register)
    refval = np.nonzero(register)[0]
    confs = list(itertools.product([0, 1], repeat=M))
    # get projection operators
    projectors = []
    for conf in confs:
        opers = [qeye(2) for i in range(N)]
        for i, c in zip(refval, list(conf)):
            opers[i] = basis(2, c).proj()
        P = tensor(opers)
        projectors += [P]
    # project
    outcomes = []
    for P in projectors:
        psip = None
        norm = 0.
        if dim2 == 1:
            psip = P*psi
            norm = psip.norm()
        elif dim2 == 2:
            psip = P*psi*P
            norm = psip.tr()
        if norm > 0. and normalize:
            psip = psip/norm
        outcomes += [psip]
    return outcomes, projectors, confs


def evolve(H, t, psi, res=200, c_ops=[]):
    opts = Options(store_final_state=True)
    times = np.linspace(0., t, res)
    result = mesolve(H, psi, times, c_ops, options=opts)
    return result.final_state
