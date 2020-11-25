#!/usr/bin/env python3
from time import time
from simpleqsimulator import QCircuit, QState
import numpy as np


# Create a 16 qubits
qc = QCircuit(1)

# Apply the two gates
qc.U(0, np.matrix([[1, 0], [0, -1j]]))

# Apply the State |0> (in the computational basis) to the circuit
qc.evaluate(QState(1, 1))

# Draw the circuit to the terminal
qc.draw()
