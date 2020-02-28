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
