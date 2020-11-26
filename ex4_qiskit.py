#!/usr/bin/env python3
import numpy as np
from qiskit import(
  QuantumCircuit,
  execute,
  Aer)
from qiskit.visualization import plot_histogram

# Use Aer's qasm_simulator
simulator = Aer.get_backend('qasm_simulator')



for gate in '1 X CNOT'.split():
    # Draw the circuit to the terminal
    print(gate)

    circuit = QuantumCircuit(2, 1)

    # Prepare the initial state |01>
    circuit.x(0)

    # Apply H to both qubits
    circuit.h(0)
    circuit.h(1)

    print('Using a %s gate as Uf' % gate)

    # Apply the gate under test
    if gate == '1':
        pass
    elif gate == 'X':
        circuit.x(0)
    elif gate == 'CNOT':
        circuit.cnot(1, 0)

    circuit.h(1)

    # Map the quantum measurement to the classical bit
    # Qubit 1 will be measured
    circuit.measure([1], [0])

    # Execute the circuit on the qasm simulator
    job = execute(circuit, simulator, shots=1000)

    # Grab results from the job
    result = job.result()

    # Print the result
    counts = result.get_counts(circuit)
    if not '0' in counts:           # result was 100% 1
        print('function f corresponding to %s is balanced' % gate)
    elif not '1' in counts:         # result was 100% 0
        print('function f corresponding to %s is constant' % gate)

    # Draw the circuit
    print(circuit.draw())
