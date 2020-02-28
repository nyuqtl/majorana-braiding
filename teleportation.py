import numpy as np

from qutip import basis, controlled_gate, tensor, expect
from qutip import sigmaz, snot, rx, ry, rz

# functions related to quantum mechanical concepts
from qm import evolve, pmeasurement
from qm import Sx, Sy, Sz

# Hamiltonian generators
from hamiltonians import constructHadamardH
from hamiltonians import constructHadamardCorr
from hamiltonians import constructCZH

# helper functions
from helpers import osum


# schedule is list of unitary operators
def scheduledUnitaryEvolution(psi0, schedule):
    psif = None
    for U in schedule:
        if psif is None:
            psif = U*psi0
        else:
            psif = U*psif
    return psif


# schedule is list of triples where first element is
# a Hamiltonian, second element is correcting unitary
# and third element is evolution time
def scheduledTimeEvolution(psi0, schedule, c_ops=[]):
    psif = None
    for H, U, t in schedule:
        if psif is None:
            psif = evolve(H, t/2., psi0, c_ops=c_ops)
        else:
            psif = evolve(H, t/2., psif, c_ops=c_ops)
        if U is not None:
            if psif.dims[1][0] == 1:
                # state vector
                psif = U*psif
            else:
                # density matrix
                psif = U.dag()*psif*U
    return psif


# takes a quantum state to be teleported
# and decay rate, produces a state of 8 qubits
# runs encoding state and teleportation
# returns a reuslting density operator
def continuousTeleportationSimulation(psi, c_ops=[]):
    N = 8
    psi0 = tensor([basis(2, 0), psi] + [basis(2, 0) for i in range(N-2)])
    schedule = []
    # encoding, Hadamards stage 1, 2
    schedule.append((
        constructHadamardH(N, range(N)),
        constructHadamardCorr(N, range(N)),
        np.pi
    ))
    # encoding, CZs, stage 3
    schedule.append((
        constructCZH(N, [0, 2, 4], [1, 3, 5]),
        None,
        np.pi
    ))
    # encoding, Hadamards, stage 4, 5
    schedule.append((
        constructHadamardH(N, [1, 3, 5]),
        constructHadamardCorr(N, [1, 3, 5]),
        np.pi
    ))
    # teleportation, stage 6
    schedule.append((
        osum([Sy(N, i) for i in [1, 2, 3, 4]]),
        None,
        np.pi/2.
    ))
    # teleportation, stage 7
    schedule.append((
        osum([Sz(N, i) for i in [1, 2, 3, 4]]),
        None,
        -np.pi/2.
    ))
    # teleportation, stage 8
    schedule.append((
        constructCZH(N, [1, 3], [2, 4]),
        None,
        np.pi
    ))
    # teleportation, stage 9
    schedule.append((
        osum([Sy(N, i) for i in [1, 2, 3, 4]]),
        None,
        -np.pi/2.
    ))
    # decoding, stage 10, 11
    schedule.append((
        constructHadamardH(N, [1, 3]),
        constructHadamardCorr(N, [1, 3]),
        np.pi
    ))
    # decoding, stage 12
    schedule.append((
        constructCZH(N, [0, 2], [1, 3]),
        None,
        np.pi
    ))
    # decoding, stage 13, 14
    schedule.append((
        constructHadamardH(N, [0, 1, 2, 3]),
        constructHadamardCorr(N, [0, 1, 2, 3]),
        np.pi
    ))
    # perform time evolution and return the final state
    return scheduledTimeEvolution(psi0, schedule, c_ops=c_ops)


def circuitTeleportationSimulation(psi, c_ops=[]):
    N = 8
    psi0 = tensor([basis(2, 0), psi] + [basis(2, 0) for i in range(N-2)])
    schedule = []
    # encoding, Hadamards stage 1, 2
    schedule += [snot(N=N, target=i) for i in range(N)]
    # encoding, CZs, stage 3
    schedule += [
        controlled_gate(sigmaz(), N=N, control=0, target=1),
        controlled_gate(sigmaz(), N=N, control=2, target=3),
        controlled_gate(sigmaz(), N=N, control=4, target=5)
    ]
    # encoding, Hadamards, stage 4, 5
    schedule += [snot(N=N, target=i) for i in [1, 3, 5]]
    # teleportation, stage 6
    schedule += [ry(np.pi/2., N=N, target=i) for i in [1, 2, 3, 4]]
    # teleportation, stage 7
    schedule += [rz(-np.pi/2., N=N, target=i) for i in [1, 2, 3, 4]]
    # teleportation, stage 8
    schedule += [
        controlled_gate(sigmaz(), N=N, control=1, target=2),
        controlled_gate(sigmaz(), N=N, control=3, target=4)
    ]
    # teleportation, stage 9
    schedule += [ry(-np.pi/2., N=N, target=i) for i in [1, 2, 3, 4]]
    # decoding, stage 10, 11
    schedule += [snot(N=N, target=i) for i in [1, 3]]
    # decoding, stage 12
    schedule += [
        controlled_gate(sigmaz(), N=N, control=0, target=1),
        controlled_gate(sigmaz(), N=N, control=2, target=3)
    ]
    # decoding, stage 13, 14
    schedule += [snot(N=N, target=i) for i in [0, 1, 2, 3]]
    # perform time evolution and return the final state
    return scheduledUnitaryEvolution(psi0, schedule)


def continuousXXBraidingCorrectionSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    for i in range(2):
        # XX braiding, stage 1
        schedule.append((
            osum([Sy(N, i) for i in [5, 6]]),
            None,
            np.pi/2.
        ))
        # XX braiding, stage 2
        schedule.append((
            osum([Sz(N, i) for i in [5, 6]]),
            None,
            -np.pi/2.
        ))
        # XX braiding, stage 3
        schedule.append((
            constructCZH(N, [5], [6]),
            None,
            np.pi
        ))
        # XX braiding, stage 4
        schedule.append((
            osum([Sy(N, i) for i in [5, 6]]),
            None,
            -np.pi/2.
        ))
    return scheduledTimeEvolution(psi0, schedule, c_ops=c_ops)


def circuitXXBraidingCorrectionSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    for i in range(2):
        # XX braiding, stage 1
        schedule += [ry(np.pi/2., N=N, target=i) for i in [5, 6]]
        # XX braiding, stage 2
        schedule += [rz(-np.pi/2., N=N, target=i) for i in [5, 6]]
        # XX braiding, stage 3
        schedule += [
            controlled_gate(sigmaz(), N=N, control=5, target=6)
        ]
        # XX braiding, stage 4
        schedule += [ry(-np.pi/2., N=N, target=i) for i in [5, 6]]
    # perform time evolution and return the final state
    return scheduledUnitaryEvolution(psi0, schedule)


def continuousZBraidingCorrectionSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    for i in range(2):
        # XX braiding, stage 1
        schedule.append((
            osum([Sx(N, i) for i in [4, 5]]),
            None,
            np.pi/2.
        ))
        # XX braiding, stage 2
        schedule.append((
            osum([Sz(N, i) for i in [4, 5]]),
            None,
            -np.pi/2.
        ))
        # XX braiding, stage 3
        schedule.append((
            constructCZH(N, [4], [5]),
            None,
            np.pi
        ))
        # XX braiding, stage 4
        schedule.append((
            osum([Sx(N, i) for i in [4, 5]]),
            None,
            -np.pi/2.
        ))
    return scheduledTimeEvolution(psi0, schedule, c_ops=c_ops)


def circuitZBraidingCorrectionSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    for i in range(2):
        # Z braiding, stage 1
        schedule += [rx(np.pi/2., N=N, target=i) for i in [4, 5]]
        # Z braiding, stage 2
        schedule += [rz(-np.pi/2., N=N, target=i) for i in [4, 5]]
        # Z braiding, stage 3
        schedule += [
            controlled_gate(sigmaz(), N=N, control=4, target=5)
        ]
        # Z braiding, stage 4
        schedule += [rx(-np.pi/2., N=N, target=i) for i in [4, 5]]
    # perform time evolution and return the final state
    return scheduledUnitaryEvolution(psi0, schedule)


def continuousDecodingSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    # decoding, stage 1, 2
    schedule.append((
        constructHadamardH(N, [5, 7]),
        constructHadamardCorr(N, [5, 7]),
        np.pi
    ))
    # decoding, stage 3
    schedule.append((
        constructCZH(N, [4, 6], [5, 7]),
        None,
        np.pi
    ))
    # decoding, stage 4, 5
    schedule.append((
        constructHadamardH(N, [4, 5, 6]),
        constructHadamardCorr(N, [4, 5, 6]),
        np.pi
    ))
    return scheduledTimeEvolution(psi0, schedule, c_ops=c_ops)


def circuitDecodingSimulation(psi0, c_ops=[]):
    N = 8
    schedule = []
    # decoding, stage 1, 2
    schedule += [snot(N=N, target=i) for i in [5, 7]]
    # decoding, stage 3
    schedule += [
        controlled_gate(sigmaz(), N=N, control=4, target=5),
        controlled_gate(sigmaz(), N=N, control=6, target=7)
    ]
    # decoding, stage 4, 5
    schedule += [snot(N=N, target=i) for i in [4, 5, 6]]
    # perform time evolution and return the final state
    return scheduledUnitaryEvolution(psi0, schedule)


def simulateTeleportation(
        psi,
        Ftel,
        FXX,
        FZ,
        Fdec,
        normalize=True, c_ops=[]):
    amZ = 0.  # post-selected amplitudes
    smZ = 0.
    amn = 0.  # all the amplitudes
    smn = 0.
    psifc = Ftel(psi, c_ops=c_ops)
    M = [True, True, True, True, False, False, False, False]
    psis, _, outs = pmeasurement(psifc, M, normalize=False)
    for psiout, out in zip(psis, outs):
        c1, c2 = out[1], out[3]
        psio = None
        if c1 == 0 and c2 == 0:
            psio = FZ(psiout, c_ops=c_ops)
        if c1 == 0 and c2 == 1:
            psio = FXX(psiout, c_ops=c_ops)
        if c1 == 1 and c2 == 0:
            psio = FXX(psiout, c_ops=c_ops)
            psio = FZ(psio, c_ops=c_ops)
        if c1 == 1 and c2 == 1:
            psio = psiout
        # apply the decoding circuit
        psif = Fdec(psio, c_ops=c_ops)
        # perform another projection on remaining qubits
        M = [False, False, False, False, True, False, True, True]
        mpsis, _, mouts = pmeasurement(psif, M, normalize=False)
        for mpsif, mout in zip(mpsis, mouts):
            # extract the teleported state
            psiexp = tensor([
                basis(2, out[0]),
                basis(2, out[1]),
                basis(2, out[2]),
                basis(2, out[3]),
                basis(2, mout[0]),
                psi,
                basis(2, mout[1]),
                basis(2, mout[2])])
            vv = None
            nrm = None
            if mpsif.dims[1][0] == 1:
                # state vector
                vv = np.abs(mpsif.overlap(psiexp))
                nrm = mpsif.norm()
            else:
                # density matrix
                vv = expect(mpsif, psiexp)
                nrm = mpsif.tr()
            if np.count_nonzero([out[0], out[2], mout[0], mout[1]]) == 0:
                amZ += vv
                smZ += nrm
            amn += vv
            smn += nrm
    ampZ = 0.
    ampn = 0.
    ampZ = amZ / smZ
    ampn = amn / smn
    return ampZ, ampn
