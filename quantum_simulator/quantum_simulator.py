#!/usr/bin/env python3
import numpy as np				# numerical py
from scipy import sparse		# scientific py


# Characters to draw the circuit
font = '┌─┐│└┘┥┝┴┬━┿■X'


class Gate:
    '''
    Object representing one gate, currently only used for drawing the circuit.
    '''

    def __init__(self, type_, controls=[], targets=[]):
        # type_: SWAP, X, Y, ...
        self.type_ = type_
        self.controls = controls
        self.targets = targets

    def draw(self, permutation, gate_width=5, gate_height=3, gate_separation_height=1):
        # Renders the gate into a buffer and returns the required space
        controls = [permutation[c] for c in self.controls]
        targets = [permutation[t] for t in self.targets]
        used_bit_positions = controls + targets
        min_qubit = min(used_bit_positions)
        max_qubit = max(used_bit_positions)
        gate_distance = gate_height + gate_separation_height
        min_row = min_qubit * gate_distance
        max_row = max_qubit * gate_distance + gate_height
        output = [''] * (max_row - min_row)

        all_rows_center = set(range(min_qubit * gate_distance - min_row + (gate_height // 2),
            max_qubit * gate_distance - min_row + (gate_height // 2)))

        def to_centers(qubits):
            # Returns the center line of the qubits
            return [qubit_number * gate_distance - min_row + (gate_height // 2) for qubit_number in qubits]

        def to_rows(qubits):
            # Returns the rows affected by the qubits
            result = []
            for qubit_number in qubits:
                for row in range(qubit_number * gate_distance, qubit_number * gate_distance + gate_height):
                    result.append(row - min_row)
            return result

        if self.type_ != 'SWAP':
            # SWAP has no label
            for control in controls:
                row_middle = control * gate_distance - min_row + (gate_height // 2)
                output[row_middle] += font[12].center(gate_width, font[10])

        # Calculate crossings
        crossings = []
        for qubit in range(min_qubit, max_qubit):
            if not qubit in controls + targets:
                row_middle = qubit * gate_distance - min_row + (gate_height // 2)
                crossings.append(qubit)
                output[row_middle] += font[11].center(gate_width, font[10])

        if self.type_ == 'SWAP':
            # Use an "X" to symbolize the Swap gate
            for target in targets:
                row_middle = target * gate_distance - min_row + (gate_height // 2)
                output[row_middle] += font[13].center(gate_width, font[10])

            # control qubit line intersects qubit
            for row in all_rows_center - set(to_centers(controls + targets + crossings)):
                output[row] += font[3].center(gate_width)

            min_row += (gate_height - 1) // 2   # Swap does not have a box and is therefore smaller
            max_row -= (gate_height - 1) // 2

        else:
            # draw target(s)
            for target in targets:
                row_start = target * gate_distance - min_row
                row_stop = gate_height + target * gate_distance - min_row - 1
                row_middle = (row_start + row_stop) // 2
                rest = set(range(row_start + 1, row_stop)) - set([row_middle])
                # draw upper part of box
                output[row_start] += font[0] + font[8 if any(c < target for c in controls) else 1].center(gate_width - 2, font[1]) + font[2]
                # draw middle part of box
                output[row_middle] += font[6] + self.type_.center(gate_width - 2) + font[7]
                # draw lower part of box
                output[row_stop] += font[4] + font[9 if any(c > target for c in controls) else 1].center(gate_width - 2, font[1]) + font[5]
                # draw rest of box
                for row in rest:
                    output[row] += font[3] + ' ' * (gate_width - 2) + font[3]

            all_rows = set(range(min_qubit * gate_distance - min_row,
                max_qubit * gate_distance - min_row + gate_height))

            for row in all_rows - set(to_rows(targets) + to_centers(controls + crossings)):
                if row in all_rows_center:
                    output[row] += font[3].center(gate_width)
                else:
                    output[row] += ' ' * gate_width

        # Returns the rendered output, the range of affected qubits and the range of affected row
        return output, set(range(min_qubit, max_qubit + 1)), set(range(min_row, max_row))


class QState:
    '''
    Class representing one 2 ** numbits-dimensional normalised vector.
    They can be multiplied with a (complex) scalar and be added and subtracted
    '''
    def __init__(self, state, numbits=None):
        self.numbits = numbits
        if isinstance(state, (list, np.ndarray)) and (numbits is None or len(state) == numbits ** 2):
            self.vector = np.array(state, dtype=np.complex)
            length = np.sqrt(np.sum(np.abs(self.vector * self.vector)))
            if np.isclose(length, 0):
                raise ValueError('Vector of length 0 not allowed')
            self.vector /= length       # normalise vector
            self.numbits = np.log2(len(self.vector))
        elif isinstance(state, (int, np.int, np.int16, np.int32, np.int64)):
            self.vector = np.zeros(2 ** numbits, dtype=np.complex)
            self.vector[state] = 1
        assert np.isclose(self.numbits % 1, 0)
        self.numbits = int(self.numbits)

    def __mul__(self, other):
        # Multiply state by a scalar
        if np.isscalar(other):
            vector = self.vector * other
            return QState(vector)
        else:
            raise ValueError('Only multiplying with a scalar is allowed')

    def __rmul__(self, other):
        # Multiply state by a scalar
        return self.__mul__(other)

    def __add__(self, other):
        # Add state "other" to state
        if isinstance(other, QState):
            return QState(self.vector + other.vector)
        else:
            raise ValueError('Only adding of QStates is alloed')

    def __radd__(self, other):
        # Add state "other" to state
        return self.__add__(other)

    def __sub__(self, other):
        # Subtract other from state
        if isinstance(other, QState):
            return QState(self.vector - other.vector)
        else:
            raise ValueError('Only subtracting of QStates is alloed')

    def __rsub__(self, other):
        # Subtract state from other
        if isinstance(other, QState):
            return QState(other.vector - self.vector)
        else:
            raise ValueError('Only subtracting of QStates is alloed')

    def __repr__(self):
        # pretty-print a state as ket
        vector_mask = ~np.isclose(self.vector, 0)
        if np.isclose(np.imag(self.vector), 0).all():
            self.vector = self.vector.real
        buffer = []
        for i, v in zip(np.arange(2 ** self.numbits)[vector_mask], self.vector[vector_mask]):
            # print resultant vector with amplitude
            buffer.append('\t %s %s' % (np.round(v, 3), '|{:0{}b}>'.format(i, self.numbits)))
        return '\n'.join(buffer)

class QCircuit:
    '''
    Class representing one quantum circuit.
    It uses the scipy sparse matrix, which is an improved version that saves only non-zero elements.
    This is only efficient for a sparsely populated matrix. However using 16 qubits, one needs a
    total of (2 ** 16) * 8 Bytes ~ 34 GByte which is far too much.
    '''

    def __init__(self, numbits=1):
        # initialiser, numbits holds the number of qubits
        assert numbits > 0
        self.numbits = numbits
        self.numstates = 2 ** self.numbits
        self.allstates = np.arange(self.numstates)

        self.gates = []
        self.gate_objects = []

    def create_single_qubit_gate(self, target, gate_matrix):
        # Add one matrix of the form
        # a b
        # c d
        assert target in range(self.numbits)

        mask_bit = 1 << target  # bit of the target
        # Get coordinates on which the single qubit matrix gets inserted
        # Coordinates where the target bit is not set
        unique = self.allstates[~np.array(np.bitwise_and(self.allstates, mask_bit), dtype=np.bool)]
        # Coordinates where the target bit was set
        unique_t = np.bitwise_or(unique, mask_bit)

        # Done in such a way that all (row, column) pairs for example (unique, unique) correspond to the
        # value gate_matrix[0, 0]

        # Rows and columns where all bits except for the target bits are equal for input and output
        rows = np.concatenate((unique, unique, unique_t, unique_t))
        cols = np.concatenate((unique, unique_t, unique, unique_t))
        #                      a,      b,        c,      d

        size4 = self.numstates // 2
        values = np.concatenate((
            np.full(size4, gate_matrix[0, 0]),
            np.full(size4, gate_matrix[0, 1]),
            np.full(size4, gate_matrix[1, 0]),
            np.full(size4, gate_matrix[1, 1])
        ))

        # Create the sparsely populated matrix
        matrix = sparse.csr_matrix((values, (rows, cols)))
        self.gates.append(matrix)

    def create_controlled_single_qubit_gate(self, control, target, gate_matrix):
        # Add one matrix of the form
        # 1 0 0 0
        # 0 1 0 0
        # 0 0 a b
        # 0 0 c d

        assert target in range(self.numbits)
        assert control in range(self.numbits)
        assert target != control

        target_bit = 1 << target  # bit of the target
        control_bit = 1 << control

        # Get coordinates on which the single qubit matrix gets inserted
        unique = np.unique(np.bitwise_and(self.allstates, self.numstates - 1 - (target_bit + control_bit)))
        unique_t = np.bitwise_or(unique, target_bit)
        unique_c = np.bitwise_or(unique, control_bit)
        unique_ct = np.bitwise_or(unique, target_bit + control_bit)

        # Done in such a way that all (row, column) pairs for example (unique, unique) correspond to the
        # value gate_matrix[0, 0]

        rows = np.concatenate((unique, unique_t, unique_c, unique_ct, unique_c, unique_ct))
        cols = np.concatenate((unique, unique_t, unique_c, unique_c, unique_ct, unique_ct))
        #                      1,      1,        a,        b,        c,         d

        size4 = self.numstates // 4
        values = np.concatenate((
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, gate_matrix[0, 0]),
            np.full(size4, gate_matrix[0, 1]),
            np.full(size4, gate_matrix[1, 0]),
            np.full(size4, gate_matrix[1, 1])
        ))

        # Create the sparsely populated matrix
        matrix = sparse.csr_matrix((values, (rows, cols)))
        self.gates.append(matrix)

    def create_double_qubit_gate(self, target1, target2, gate_matrix):
        # Add one matrix of the form
        # a b c d
        # e f g h
        # i j k l
        # m n o p

        assert target1 in range(self.numbits)
        assert target2 in range(self.numbits)
        assert target1 != target2

        target1_bit = 1 << target1  # bit of the target
        target2_bit = 1 << target2

        # Get coordinates on which the single qubit matrix gets inserted
        unique = np.unique(np.bitwise_and(self.allstates, self.numstates - 1 - (target1_bit + target2_bit)))
        unique_t1 = np.bitwise_or(unique, target1_bit)
        unique_t2 = np.bitwise_or(unique, target2_bit)
        unique_t1t2 = np.bitwise_or(unique, target1_bit + target2_bit)

        # Done in such a way that all (row, column) pairs for example (unique, unique) correspond to the
        # value gate_matrix[0, 0]

        rows = np.concatenate((
            unique, unique, unique, unique,
            unique_t1, unique_t1, unique_t1, unique_t1,
            unique_t2, unique_t2, unique_t2, unique_t2,
            unique_t1t2, unique_t1t2, unique_t1t2, unique_t1t2
        ))
        cols = np.concatenate((
            unique, unique_t1, unique_t2, unique_t1t2,
            unique, unique_t1, unique_t2, unique_t1t2,
            unique, unique_t1, unique_t2, unique_t1t2,
            unique, unique_t1, unique_t2, unique_t1t2,
        ))

        size4 = self.numstates // 4
        values = np.concatenate((
            np.full(size4, gate_matrix[0, 0]),
            np.full(size4, gate_matrix[0, 1]),
            np.full(size4, gate_matrix[0, 2]),
            np.full(size4, gate_matrix[0, 3]),
            np.full(size4, gate_matrix[1, 0]),
            np.full(size4, gate_matrix[1, 1]),
            np.full(size4, gate_matrix[1, 2]),
            np.full(size4, gate_matrix[1, 3]),
            np.full(size4, gate_matrix[2, 0]),
            np.full(size4, gate_matrix[2, 1]),
            np.full(size4, gate_matrix[2, 2]),
            np.full(size4, gate_matrix[2, 3]),
            np.full(size4, gate_matrix[3, 0]),
            np.full(size4, gate_matrix[3, 1]),
            np.full(size4, gate_matrix[3, 2]),
            np.full(size4, gate_matrix[3, 3])
        ))

        # Create the sparsely populated matrix
        matrix = sparse.csr_matrix((values, (rows, cols)))
        self.gates.append(matrix)

    def create_double_controlled_single_qubit_gate(self, control1, control2, target, gate_matrix):
        # Add one matrix of the form
        # 1 0 0 0 0 0 0 0
        # 0 1 0 0 0 0 0 0
        # 0 0 1 0 0 0 0 0
        # 0 0 0 1 0 0 0 0
        # 0 0 0 0 1 0 0 0
        # 0 0 0 0 0 1 0 0
        # 0 0 0 0 0 0 a b
        # 0 0 0 0 0 0 c d

        assert target in range(self.numbits)
        assert control1 in range(self.numbits)
        assert control2 in range(self.numbits)
        assert target != control1
        assert target != control2
        assert control1 != control2

        target_bit = 1 << target  # bit of the target
        control1_bit = 1 << control1
        control2_bit = 1 << control2

        # Get coordinates on which the single qubit matrix gets inserted
        unique = np.unique(np.bitwise_and(self.allstates, self.numstates - 1 - (target_bit + control1_bit + control2_bit)))
        unique_t = np.bitwise_or(unique, target_bit)
        unique_c1 = np.bitwise_or(unique, control1_bit)
        unique_c1t = np.bitwise_or(unique, target_bit + control1_bit)
        unique_c2 = np.bitwise_or(unique, control2_bit)
        unique_c2t = np.bitwise_or(unique, target_bit + control2_bit)
        unique_c1c2 = np.bitwise_or(unique, control1_bit + control2_bit)
        unique_c1c2t = np.bitwise_or(unique, target_bit + control2_bit + control1_bit)

        # Done in such a way that all (row, column) pairs for example (unique, unique) correspond to the
        # value gate_matrix[0, 0]

        rows = np.concatenate((unique, unique_t, unique_c1, unique_c1t, unique_c2, unique_c2t, unique_c1c2, unique_c1c2t, unique_c1c2, unique_c1c2t))
        cols = np.concatenate((unique, unique_t, unique_c1, unique_c1t, unique_c2, unique_c2t, unique_c1c2, unique_c1c2, unique_c1c2t, unique_c1c2t))

        size4 = self.numstates // 8
        values = np.concatenate((
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, 1),
            np.full(size4, gate_matrix[0, 0]),
            np.full(size4, gate_matrix[0, 1]),
            np.full(size4, gate_matrix[1, 0]),
            np.full(size4, gate_matrix[1, 1])
        ))

        # Create the sparsely populated matrix
        matrix = sparse.csr_matrix((values, (rows, cols)))
        self.gates.append(matrix)

    def identity(self):
        # Identity gate, used as a placeholder
        self.create_single_qubit_gate(0, np.array([[1, 0], [0, 1]]))
        self.gate_objects.append(Gate('identity'))

    def H(self, target):
        # Hadamard gate
        self.create_single_qubit_gate(target, np.array([[1, 1], [1, -1]]) / np.sqrt(2))
        self.gate_objects.append(Gate('H', targets=[target]))

    def X(self, target):
        # Pauli-X gate
        self.create_single_qubit_gate(target, np.array([[0, 1], [1, 0]]))
        self.gate_objects.append(Gate('X', targets=[target]))

    def Y(self, target):
        # Pauli-Y gate
        self.create_single_qubit_gate(target, np.array([[0, -1j], [1j, 0]]))
        self.gate_objects.append(Gate('Y', targets=[target]))

    def CY(self, control, target):
        # Pauli-Y gate
        self.create_controlled_single_qubit_gate(control, target, np.array([[0, -1j], [1j, 0]]))
        self.gate_objects.append(Gate('Y', targets=[target], controls=[control]))

    def Z(self, target):
        # Pauli-Z gate
        self.create_single_qubit_gate(target, np.array([[1, 0], [0, -1]]))
        self.gate_objects.append(Gate('Z', targets=[target]))

    def S(self, target):
        # Phase gate
        self.create_single_qubit_gate(target, np.array([[1, 0], [0, 1j]]))
        self.gate_objects.append(Gate('S', targets=[target]))

    def R(self, angle, target):
        # Rotation by exp(i angle)
        self.create_single_qubit_gate(target, np.array([[1, 0], [0, np.exp(1j * angle)]]))
        self.gate_objects.append(Gate('R', targets=[target]))

    def T(self, target, invert=False):
        # PI / 8 gate
        angle = -np.pi / 4 if invert else np.pi / 4
        self.create_single_qubit_gate(target, np.array([[1, 0], [0, np.exp(1j * angle)]]))
        self.gate_objects.append(Gate('T', targets=[target]))

    def U(self, target, matrix, name='U'):
        # arbitrary unitary gate
        assert isinstance(matrix, np.matrix)            # only np.matrix and derived classes allowed
        assert matrix.shape == (2, 2)                   # only 2x2 matricies are allowed
        assert np.allclose(np.eye(matrix.shape[0]), matrix.H * matrix), 'matrix not unitary'  # matrix has to be unitary
        self.create_single_qubit_gate(target, matrix)
        self.gate_objects.append(Gate(name, targets=[target]))

    def CR(self, k, control, target, invert=False):
        # Controlled rotation by exp(2 pi i / 2 ^ k)
        angle = -2 * np.pi / 2 ** k if invert else 2 * np.pi / 2 ** k
        self.create_controlled_single_qubit_gate(control, target, np.array([[1, 0], [0, np.exp(1j * angle)]]))
        self.gate_objects.append(Gate('R%d%s' % (k, '†' if invert else ''), controls=[control], targets=[target]))

    def cnot(self, control, target):
        # PI / 8 gate
        self.create_controlled_single_qubit_gate(control, target, np.array([[0, 1], [1, 0]]))
        self.gate_objects.append(Gate('X', controls=[control], targets=[target]))

    def swap(self, target1, target2):
        # swap gate
        self.create_double_qubit_gate(target1, target2, np.array([[1, 0, 0, 0], [0, 0, 1, 0,], [0, 1, 0, 0], [0, 0, 0, 1]]))
        self.gate_objects.append(Gate('SWAP', targets=[target1, target2]))

    def toffoli(self, control1, control2, target):
        # toffoli gate
        self.create_double_controlled_single_qubit_gate(control1, control2, target, np.array([[0, 1], [1, 0]]))
        self.gate_objects.append(Gate('X', controls=[control1, control2], targets=[target]))

    def compile_matrix(self):
        # multiplies matricies of the gates so that it has only to be done once
        if len(self.gates) == 0:
            self.identity()
        m = self.gates[-1]
        for matrix in self.gates[::-1][1::]:
            m = m.__matmul__(matrix)
        return m

    def evaluate(self, state=None):
        # Evaluate all (None), one (int) or a list of states (list)

        # Handle different argument types
        if state == None:
            if self.numbits > 6:
                raise ValueError('Please don\'t do this!')
            states = [QState(state, numbits=self.numbits) for state in self.allstates]
        elif isinstance(state, int):
            assert 0 <= state < self.numstates
            states = [QState(state, numbits=self.numbits)]
        elif isinstance(state, QState):
            states = [state]
        elif isinstance(state, list):
            states = state
        else:
            raise ValueError('State of type %s not compatible' % type(state))

        # create matrix
        #matrix = self.compile_matrix()

        # evaluate all states
        for state in states:
            # create input vector

            print('input:')
            print(state)

            # create output vector
            for matrix in self.gates:
                state = QState(matrix.dot(state.vector))

            print('output:')
            print(state)

    def pprint(self):
        # pretty-print the matrix, if it is reasonable
        if self.numbits > 6:
            raise ValueError('Please don\'t do this!')
        matrix = self.compile_matrix().toarray()
        if np.isclose(np.imag(matrix), 0).all():
            matrix = matrix.real
        print(np.round(matrix, 3))

    def qubit_to_row_range(self, numbers, gate_width=None, gate_height=3, gate_separation_height=1):
        gate_distance = gate_height + gate_separation_height
        if len(numbers) == 0:
            return []
        return list(range(min(numbers) * gate_distance, max(numbers) * gate_distance + gate_height))

    def draw(self, permutation=None):
        '''
        Draws the circuit to the terminal
        '''
        gate_width = 7
        gate_height = 3
        gate_separation_width = 2
        gate_separation_height = 1
        gate_distance = gate_height + gate_separation_height

        assert gate_height % 2 == 1

        params = {
            'gate_width': gate_width,
            'gate_height': gate_height,
            'gate_separation_height': gate_separation_height,
        }

        height = gate_height * self.numbits + gate_separation_height * (self.numbits - 1)

        if permutation == None:
            permutation = {val: val for val in range(self.numbits)}

        separator = [' ' * gate_separation_width] * height
        label_width = len(str(self.numbits)) + 1
        labels = [' ' * label_width] * height
        for qubit in range(self.numbits):
            row = qubit * gate_distance + (gate_height - 1) // 2
            separator[row] = font[10] * gate_separation_width
            labels[row] = str(permutation[qubit]).ljust(label_width)

        output = labels[:]
        slice_output = separator[:]
        slice_qubits = set()
        slice_rows = set()
        gate_qubits = set()
        gate_rows = set()
        all_rows = set()

        # extra None to force render last slice
        for gate in self.gate_objects + [None]:
            create_new_slice = True
            if gate != None:
                # Render gate
                gate_output, gate_qubits, gate_rows = gate.draw(permutation, **params)

                # check if it would fit into current slice
                create_new_slice = len(slice_rows & gate_rows) > 0

            if create_new_slice: # at least one gate position does not fit
                # Add wires to empty gate positions
                for qubit_num in set(range(self.numbits)) - slice_qubits:
                    row = qubit_num * gate_distance + (gate_height - 1) // 2
                    slice_output[row] += font[10] * gate_width
                    slice_rows.add(row)     # Mark row as filled

                # Fill empty rows with blanks
                for row in set(range(height)) - slice_rows:
                    slice_output[row] += ' ' * gate_width

                # Write slice to the output
                for row, row_output in enumerate(slice_output):
                    output[row] += row_output

                all_rows |= slice_rows

                # Create new slice
                slice_output = separator[:]
                slice_qubits = set()
                slice_rows = set()

            # Now the gate should fit
            for gate_row, row in enumerate(self.qubit_to_row_range(gate_qubits, **params)):
                slice_output[row] += gate_output[gate_row]

            slice_qubits |= gate_qubits
            slice_rows |= gate_rows

        # terminate with short wires for symmetry
        for row, char in enumerate(separator):
            output[row] += char

        # Remove empty lines
        empty = set(range(height)) - all_rows
        gsh2 = gate_separation_height // 2
        remove_lines = set(l for l in empty if set(range(l - gsh2, l + (gate_separation_height - gsh2) + 1)).issubset(empty))
        print('\n'.join(o for i, o in enumerate(output) if i not in remove_lines))


# USAGE
# create a quantum circuit object with n qubits, e.g.
# qc = QCircuit(3)

# apply some gates
# qc.H(0)               # Hadamard on the 0-th gate
# qc.cnot(0, 1)         # cnot with control 0 and target 1

# draw the circuit to the terminal
# qc.draw()

# Print the matrix
# qc.pprint()

# Evaluate the matrix for some states, all 3 ^ 2 states
# qc.evaluate()

# Only evaluate states 000, 001, 010
# qc.evaluate([QState(i, 3) for i in range(3)])
