from time import time
from quantum_simulator.quantum_simulator import QC

t = time()
qc = QC(16)
qc.X(1)
qc.H(14)
print(time() - t)
qc.evaluate(0)
qc.draw()
