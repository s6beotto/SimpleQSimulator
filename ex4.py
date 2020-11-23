#!/usr/bin/env python3
from quantum_simulator.quantum_simulator import QCircuit, QState
import numpy as np


for gate in '1 X CNOT XCNOT'.split():
    # Create a 2 qubit-circuit
    qc = QCircuit(2)

    qc.H(0)
    qc.H(1)

    # Apply the gate under test
    if gate == '1':
        # Just for clarity in the drawn circuit
        qc.U(0, np.matrix([[1, 0], [0, 1]]), name='1')
    elif gate == 'X':
        qc.X(0)
    elif gate == 'CNOT':
        qc.cnot(1, 0)

    qc.H(1)

    # Apply the State |01> (in the computational basis) to the circuit
    qc.evaluate(QState(1, 2))
    value, state = qc.measure(QState(1, 2), bits=1)
    print('measured:', value)
    print('function %s is %s' % (gate, 'balanced' if value else 'constant'))

    # Draw the circuit to the terminal
    qc.draw()
