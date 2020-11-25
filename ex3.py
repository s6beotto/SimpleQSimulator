#!/usr/bin/env python3
from time import time
from simpleqsimulator import QCircuit, QState


# Save start time to calculate execution speed
t = time()

# Create a 16 qubits
qc = QCircuit(16)

# Apply the two gates
qc.X(1)
qc.H(14)

# print the execution time
print(time() - t)

# Apply the State |0> (in the computational basis) to the circuit
qc.evaluate(QState(0, 16))

# Draw the circuit to the terminal
qc.draw()
