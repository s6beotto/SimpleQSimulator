from IPython import embed
from simpleqsimulator import QState, QCircuit


def start():
    print('Loading IPython with imported library')
    print('Create a circuit as follows:')
    print('qc = QCircuit(2)')
    print('Apply some gates:')
    print('qc.H(0)')
    print('qc.cnot(0, 1)')
    print('Draw it:')
    print('qc.draw()')
    embed(colors="neutral")
