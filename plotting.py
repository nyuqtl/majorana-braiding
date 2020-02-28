import numpy as np
import argparse

import matplotlib.pyplot as plt

description = '\n'.join([
    'Majorana braiding circuit simulation',
    'using Lindblad dynamics with decoherence,',
    'by Marek Narozniak (c) GPL-3.0'
])
parser = argparse.ArgumentParser(description=description)
parser.add_argument('input', type=str, help='generated numerical results file')
parser.add_argument('output', type=str, help='output plot file')

args = parser.parse_args()

results = np.load(args.input)

gs = results[0, :]

f0000s = results[1, :]
falls = results[2, :]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('')
ax.set_ylabel('fidelity')
ax.set_xlabel('$\\gamma$')
ax.plot(gs, f0000s, label='post-selected')
ax.plot(gs, falls, label='general')
ax.legend()
ax.grid()
fig.savefig(args.output)
