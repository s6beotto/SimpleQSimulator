from time import time
from quantum_simulator.quantum_simulator import QC, QState

t = time()
qc = QC(16)
qc.X(1)
qc.H(14)
print(time() - t)
qc.evaluate(QState(0, 16))
qc.draw()
