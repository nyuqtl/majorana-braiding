[![Build Status](https://travis-ci.com/nyuqtl/majorana-teleportation.svg?branch=master)](https://travis-ci.com/nyuqtl/majorana-teleportation)

# Teleporting Majorana Zero Modes

We derived a quantum circuit which "emulates" the quantum teleportation in 1D topological regime in which quantum information is represented by Majorana Zero Modes.

This repository contains numerical simulation of this circuit with and without decoherence (using Lindblad dynamics and unitary evolution respectfully).

## Usage

There are multiple features of this repository, most important is testing the correctness of our framework but also generating fidelity plots.

### Running tests

To run tests locally simply execute the following commands

```sh
python -m pytest braids_tests.py
```

Mind the fact that full tests will take about 20 minutes as our system consists of eight qubits.

### Generating plots

Execute following command and follow the instructions

```
python fidelity.py --help
```

This command generates the datapoints for the plot, which itself can be generated using following command

```
python plotting.py
```

Again, run `--help` to check the usage.

## Results

Here are some better resolution fidelity plots allowing to compare fidelity with and without the error detection mechanism.

### Teleporting |0> state

![psi=|0>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_z0.png?raw=true "|0>")

### Teleporting |1> state

![psi=|1>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_z1.png?raw=true "|1>")

### Teleporting |+> state

![psi=|+>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_xp.png?raw=true "|+>")

### Teleporting |-> state

![psi=|->](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_xm.png?raw=true "|->")

### Teleporting |+i> state

![psi=|+i>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_yp.png?raw=true "|+i>")

### Teleporting |-i> state

![psi=|-i>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_ym.png?raw=true "|-i>")

### Teleporting a random pure state

![psi=|rnd>](https://github.com/nyuqtl/majorana-teleportation/blob/master/plots/res_hpc_rnd.png?raw=true "|rnd>")
