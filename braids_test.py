import numpy as np

from qutip import basis, snot, controlled_gate, sigmaz, tensor, rand_ket

# functions related to quantum mechanical concepts
from qm import bloch
from qm import deriveUnitary
from qm import evolve
from qm import Sz

# Hamiltonian generators
from hamiltonians import constructHadamardH
from hamiltonians import constructHadamardCorr
from hamiltonians import constructCZH

# helper functions
from helpers import are_close

# functions that generate our circuit and time evolutions
from teleportation import continuousTeleportationSimulation
from teleportation import continuousXXBraidingCorrectionSimulation
from teleportation import continuousZBraidingCorrectionSimulation
from teleportation import continuousDecodingSimulation
from teleportation import circuitTeleportationSimulation
from teleportation import circuitXXBraidingCorrectionSimulation
from teleportation import circuitZBraidingCorrectionSimulation
from teleportation import circuitDecodingSimulation
from teleportation import simulateTeleportation

z0, z1, xp, xm, yp, ym = bloch(basis(2, 0), basis(2, 1))


class TestCZandH(object):
    def testHadamardUnitary(self):
        H = constructHadamardH(1, [0])
        U = snot()
        U2 = deriveUnitary(H, np.pi)
        Ucorr = constructHadamardCorr(1, [0])
        assert Ucorr*U2 == U
        assert H.isherm

    def testHadamardUnitaryMulti(self):
        H = constructHadamardH(2, [0, 1])
        UH1 = snot(N=2, target=0)
        UH2 = snot(N=2, target=1)
        U2 = deriveUnitary(H, np.pi)
        Ucorr = constructHadamardCorr(2, [0, 1])
        assert Ucorr*U2 == UH1*UH2
        assert H.isherm

    def testHadamardHamiltonian(self):
        t = np.pi
        H = constructHadamardH(1, [0])
        U = snot()
        Ucorr = constructHadamardCorr(1, [0])
        for psi0 in [z0, z1, xp, xm, yp, ym, rand_ket(2)]:
            psiu = U*psi0
            psif = evolve(H, t/2., psi0)
            psic = Ucorr*psif
            overl = psiu.overlap(psic)
            assert are_close(overl, 1.)
        assert H.isherm

    def testCZUnitary(self):
        H = constructCZH(2, [0], [1])
        U = controlled_gate(sigmaz(), N=2, control=0, target=1)
        U2 = deriveUnitary(H, np.pi)
        assert U == U2
        assert H.isherm

    def testCZHamiltonian(self):
        t = np.pi
        H = constructCZH(2, [0], [1])
        U = controlled_gate(sigmaz(), N=2, control=0, target=1)
        states = [z0, z1, xp, xm, yp, ym, rand_ket(2)]
        for psiA in states:
            for psiB in states:
                psi0 = tensor([psiA, psiB])
                psiu = U*psi0
                psic = evolve(H, t/2., psi0)
                overl = psiu.overlap(psic)
                assert are_close(overl, 1., atol=1e-04)
        assert H.isherm

    def testCZHamiltonianMore(self):
        t = np.pi
        H1 = constructCZH(4, [0], [1])
        H2 = constructCZH(4, [2], [3])
        U1 = controlled_gate(sigmaz(), N=4, control=0, target=1)
        U2 = controlled_gate(sigmaz(), N=4, control=2, target=3)
        for psiA in [z0, z1]:
            for psiB in [xp, xm]:
                for psiC in [z0, z1]:
                    for psiD in [xp, xm]:
                        psi0 = tensor([psiA, psiB, psiC, psiD])
                        psiu = U2*U1*psi0
                        psic = evolve(H1, t/2., psi0)
                        psic = evolve(H2, t/2., psic)
                        overl = psiu.overlap(psic)
                        assert are_close(overl, 1., atol=1e-04)
        assert H1.isherm
        assert H2.isherm


class TestCompareCircuitEvolutions(object):
    def test_teleportation_component(self):
        states = [z0, z1, xp, xm, yp, ym, rand_ket(2)]
        for psi in states:
            psift = continuousTeleportationSimulation(psi, 0.)
            psifc = circuitTeleportationSimulation(psi)
            overl = psifc.overlap(psift)
            assert are_close(overl, 1., atol=1e-04)

    def test_unitary_teleportation(self):
        states = [z0, z1, xp, xm, yp, ym, rand_ket(2)]
        for psi in states:
            fidelity0000, fidelity = simulateTeleportation(
                psi,
                circuitTeleportationSimulation,
                circuitXXBraidingCorrectionSimulation,
                circuitZBraidingCorrectionSimulation,
                circuitDecodingSimulation)
            assert are_close(fidelity0000, 1.)
            assert are_close(fidelity, 1.)

    def test_time_teleportation(self):
        states = [z0, z1, xp, xm, yp, ym, rand_ket(2)]
        for psi in states:
            fidelity0000, fidelity = simulateTeleportation(
                psi,
                continuousTeleportationSimulation,
                continuousXXBraidingCorrectionSimulation,
                continuousZBraidingCorrectionSimulation,
                continuousDecodingSimulation)
            assert are_close(fidelity0000, 1.)
            assert fidelity < fidelity0000

    def test_noisy_time_teleportation(self):
        N = 8
        gamma = 0.05
        c_ops = [np.sqrt(gamma)*Sz(N, i) for i in range(N)]
        states = [z0, z1, xp, xm, yp, ym, rand_ket(2)]
        for psi in states:
            fidelity0000, fidelity = simulateTeleportation(
                psi,
                continuousTeleportationSimulation,
                continuousXXBraidingCorrectionSimulation,
                continuousZBraidingCorrectionSimulation,
                continuousDecodingSimulation,
                c_ops=c_ops)
            assert fidelity0000 > fidelity
