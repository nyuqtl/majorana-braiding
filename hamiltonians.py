import numpy as np

from qutip import qeye, sigmaz, tensor

from qm import Is, Sx, Sz, deriveUnitary

from helpers import osum


def constructHadamardH(N, targets):
    return osum([(Sx(N, i)+Sz(N, i))/np.sqrt(2.) for i in targets])


def constructHadamardCorr(N, targets):
    II = Is(N)
    U = deriveUnitary(qeye(2), -np.pi)
    for i in targets:
        II[i] = U
    return tensor(II)


def constructCZH(N, controls, targets):
    tot = None
    HC = 0.5*(qeye(2)-sigmaz())
    HT = sigmaz()-qeye(2)
    for control, target in zip(controls, targets):
        II = Is(N)
        II[control] = HC
        II[target] = HT
        term = tensor(II)
        if tot is None:
            tot = term
        else:
            tot += term
    return tot
