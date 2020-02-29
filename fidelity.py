import numpy as np
import argparse

from tqdm import tqdm
from qutip import basis, rand_ket

# functions related to quantum mechanical concepts
from qm import bloch
from qm import Sz

# functions that generate our circuit and time evolutions
from teleportation import continuousTeleportationSimulation
from teleportation import continuousXXBraidingCorrectionSimulation
from teleportation import continuousZBraidingCorrectionSimulation
from teleportation import continuousDecodingSimulation
from teleportation import simulateTeleportation

description = '\n'.join([
    'Majorana braiding circuit simulation',
    'using Lindblad dynamics with decoherence,',
    'by Marek Narozniak (c) GPL-3.0',
    '',
    'For |psi> argument use z0, z1, xp, xm, yp, ym, rnd strings'
])
parser = argparse.ArgumentParser(description=description)
parser.add_argument('psi', type=str, help='state to be teleported')
parser.add_argument('output', type=str, help='path to output file')
parser.add_argument('--res', type=int, default=3, help='number of decoherence runs')
parser.add_argument('--gamma', type=float, default=1.0, help='Maximum decay rate')

args = parser.parse_args()

if args.psi not in ['z0', 'z1', 'xp', 'xm', 'yp', 'ym', 'rnd']:
    raise Exception('Your input state to be teleported must be one of: z0, z1, xp, xm, yp, ym, rnd')

z0, z1, xp, xm, yp, ym = bloch(basis(2, 0), basis(2, 1))
inp = {
    'z0': z0,
    'z1': z1,
    'xp': xp,
    'xm': xm,
    'yp': yp,
    'ym': ym,
    'rnd': rand_ket(2)
}


N = 8
nb = args.res
gmax = args.gamma
gs = np.linspace(0., gmax, nb)

psi = inp[args.psi]

results = np.zeros((3, nb))
results[0, :] = gs

for i, gamma in enumerate(tqdm(gs)):
    c_ops = [np.sqrt(gamma)*Sz(N, j) for j in range(N)]
    f0000, fall = simulateTeleportation(
        psi,
        continuousTeleportationSimulation,
        continuousXXBraidingCorrectionSimulation,
        continuousZBraidingCorrectionSimulation,
        continuousDecodingSimulation,
        c_ops=c_ops)
    results[1, i] = f0000
    results[2, i] = fall

np.save(args.output, np.array(results))
