# Simple Quantum Simulator

## Classes
There are 2 classes in this project: `QState` and `QCircuit`.
`QCircuit` implements a quantum circuit containing a number of gates that can be applied.
The argument gives the number of qubits. These objects can be printed as ASCII-Art (although using unicode) to the terminal and the complete matrix can be calculated.
`QState` represents a normalised state vector, which can be added together and multiplied by scalars.
A state `QState(0, 4)` returns the pure state 0 in a 4 qubit basis.
As usual indices are starting with 0.

## Applying gates
Gates are added to the QCircuit by calling the corresponding function with the (control and) target qubit(s) as parameters. Some gates can be inverted, as the CR-gate. An arbitrary matrix can be applied using the `U` function.


## Usage

```python
#!/usr/bin/env python3
from time import time
from quantum_simulator.quantum_simulator import QCircuit, QState


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
```

which returns
```
0.01157236099243164
input:
         1.0 |0000000000000000>
output:
         0.707 |0000000000000010>
         0.707 |0100000000000010>
              
0  ━━━━━━━━━━━
              
     ┌─────┐  
1  ━━┥  X  ┝━━
     └─────┘  
              
2  ━━━━━━━━━━━
              
3  ━━━━━━━━━━━
              
4  ━━━━━━━━━━━
              
5  ━━━━━━━━━━━
              
6  ━━━━━━━━━━━
              
7  ━━━━━━━━━━━
              
8  ━━━━━━━━━━━
              
9  ━━━━━━━━━━━
              
10 ━━━━━━━━━━━
              
11 ━━━━━━━━━━━
              
12 ━━━━━━━━━━━
              
13 ━━━━━━━━━━━
              
     ┌─────┐  
14 ━━┥  H  ┝━━
     └─────┘  
              
15 ━━━━━━━━━━━
              
```
