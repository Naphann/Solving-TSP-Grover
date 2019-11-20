import qiskit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import Aer
from qiskit.circuit import Qubit

import numpy as np
from math import floor, ceil

def RTL(qc, a, b, c):
    ## fig 3 dashed
    qc.rccx(a, b, c)
    
def RTL_inv(qc, a, b, c):
    RTL(qc, a, b, c)

def RTS(qc, a, b, c):
    ## fig 3 gates 2-6
    qc.h(c)
    qc.t(c)
    qc.cx(b, c)
    qc.tdg(c)
    qc.cx(a, c)
    
def RTS_inv(qc, a, b, c):
    qc.cx(a, c)
    qc.t(c)
    qc.cx(b, c)
    qc.tdg(c)
    qc.h(c)

def SRTS(qc, a, b, c):
    ## circuit 3 dashed
    qc.h(c)
    qc.cx(c, b)
    qc.tdg(b)
    qc.cx(a, b)
    qc.t(b)
    qc.cx(c, b)
    qc.tdg(b)
    qc.cx(a, b)
    qc.t(b)
    
def SRTS_inv(qc, a, b, c):
    qc.tdg(b)
    qc.cx(a, b)
    qc.t(b)
    qc.cx(c, b)
    qc.tdg(b)
    qc.cx(a, b)
    qc.t(b)
    qc.cx(c, b)
    qc.h(c)

def RT4L(qc, a, b, c, d):
    ## fig 4
    qc.rcccx(a, b, c, d)
    
def RT4L_inv(qc, a, b, c, d):
    qc.h(d)
    qc.t(d)
    qc.cx(c, d)
    qc.tdg(d)
    qc.h(d)
    qc.t(d)
    qc.cx(b, d)
    qc.tdg(d)
    qc.cx(a, d)
    qc.t(d)
    qc.cx(b, d)
    qc.tdg(d)
    qc.cx(a, d)
    qc.h(d)
    qc.t(d)
    qc.cx(c, d)
    qc.tdg(d)
    qc.h(d)

def RT4S(qc, a, b, c, d):
    ## fig 4 dashed
    qc.h(d)
    qc.t(d)
    qc.cx(c, d)
    qc.tdg(d)
    qc.h(d)
    qc.cx(a, d)
    qc.t(d)
    qc.cx(b, d)
    qc.tdg(d)
    qc.cx(a, d)
    
def RT4S_inv(qc, a, b, c, d):
    qc.cx(a, d)
    qc.t(d)
    qc.cx(b, d)
    qc.tdg(d)
    qc.cx(a, d)
    qc.h(d)
    qc.t(d)
    qc.cx(c, d)
    qc.tdg(d)
    qc.h(d)

def apply_mct_clean(self, controls, target, ancilla_register):
    if len(controls) < 3:
        raise ValueError("there's something wrong")
    
    n = len(controls)
    ancilla = ancilla_register[:ceil((n-2)/2)]
        
    if n == 3:
        # TODO: Check for ancilla length
        self.rccx(controls[0], controls[1], ancilla[0])
        self.ccx(controls[2], ancilla[0], target)
        self.rccx(controls[0], controls[1], ancilla[0])
        return
    
    if n == 4:
        # TODO: Check for ancilla length
        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
        self.ccx(controls[3], ancilla[0], target)
        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
        return
    
    ## when controls >= 5
    
    if n % 2 == 0:
        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
        self.barrier()
        anc_idx = 1

        for i in range(3, n-1, 2):
    #         print('i = /{}'.format(i))
            self.rcccx(controls[i], controls[i+1], ancilla[anc_idx-1], ancilla[anc_idx])
            self.barrier()
            anc_idx += 1
        if (n-3)%2 == 1:
            self.ccx(controls[-1], ancilla[-1], target)
            self.barrier()
        else:
            self.rccx(controls[-2], ancilla[-2], ancilla[-1])
            self.barrier()
            self.ccx(controls[-1], ancilla[-1], target)
            self.barrier()
            self.rccx(controls[-2], ancilla[-2], ancilla[-1])
            self.barrier()
        for i in reversed(range(3, n-1, 2)):
            anc_idx -= 1
            self.rcccx(controls[i], controls[i+1], ancilla[anc_idx-1], ancilla[anc_idx])
            self.barrier()

        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
    else:
        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
        self.barrier()
        anc_idx = 1

        for i in range(3, n-3, 2):
    #         print('i = /{}'.format(i))
            self.rcccx(controls[i], controls[i+1], ancilla[anc_idx-1], ancilla[anc_idx])
            self.barrier()
            anc_idx += 1
        if (n-3)%2 == 1:
            self.ccx(controls[-1], ancilla[-1], target)
            self.barrier()
        else:
            self.rccx(controls[-2], ancilla[-2], ancilla[-1])
            self.barrier()
            self.ccx(controls[-1], ancilla[-1], target)
            self.barrier()
            self.rccx(controls[-2], ancilla[-2], ancilla[-1])
            self.barrier()
        for i in reversed(range(3, n-3, 2)):
            anc_idx -= 1
            self.rcccx(controls[i], controls[i+1], ancilla[anc_idx-1], ancilla[anc_idx])
            self.barrier()

        self.rcccx(controls[0], controls[1], controls[2], ancilla[0])
    
def apply_mct_dirty(self, controls, target, ancilla):
    # TODO: check controls to be list of bits or register
    if len(controls) == 1:
        self.cx(controls[0], target)
        return
    if len(controls) == 2:
        self.ccx(controls[0], controls[1], target)
        return
    
    if len(controls) == 3:
        SRTS(self, controls[2], ancilla[0], target)
        RTL(self, controls[0], controls[1], ancilla[0])
        SRTS_inv(self, controls[2], ancilla[0], target)
        RTL_inv(self, controls[0], controls[1], ancilla[0])
        return

    n = len(controls)
    anc = ancilla[:ceil((n-2)/2)]
    
    
    SRTS(self, controls[-1], anc[-1], target)
    qc.barrier()
    
    if (n-4)%2 == 0:
        a_idx = 1
        for i in reversed(range(floor((n-4)/2))):
            RT4S(self, anc[a_idx - 1 + i], controls[2*i+3], controls[2*i+4], anc[a_idx + i])
            qc.barrier()
    else:
        a_idx = 2
        for i in reversed(range(floor((n-4)/2))):
            RT4S(self, anc[a_idx - 1 + i], controls[2*i+4], controls[2*i+5], anc[a_idx + i])
            qc.barrier()
        RTS(self, anc[0], controls[3], anc[1])
        qc.barrier()
    
    RT4L(self, controls[0], controls[1], controls[2], anc[0])
    qc.barrier()
    
    if (n-4)%2 == 0:
        a_idx = 1
        for i in (range(floor((n-4)/2))):
            RT4S_inv(self, anc[a_idx - 1 + i], controls[2*i+3], controls[2*i+4], anc[a_idx + i])
            qc.barrier()
    else:
        a_idx = 2
        RTS_inv(self, anc[0], controls[3], anc[1])
        qc.barrier()
        for i in (range(floor((n-4)/2))):
            RT4S_inv(self, anc[a_idx - 1 + i], controls[2*i+4], controls[2*i+5], anc[a_idx + i])
            qc.barrier()
            
    SRTS_inv(self, controls[-1], anc[-1], target)
    qc.barrier()
    
    ## SAME AS ABOVE
    if (n-4)%2 == 0:
        a_idx = 1
        for i in reversed(range(floor((n-4)/2))):
            RT4S(self, anc[a_idx - 1 + i], controls[2*i+3], controls[2*i+4], anc[a_idx + i])
            qc.barrier()
    else:
        a_idx = 2
        for i in reversed(range(floor((n-4)/2))):
            RT4S(self, anc[a_idx - 1 + i], controls[2*i+4], controls[2*i+5], anc[a_idx + i])
            qc.barrier()
        RTS(self, anc[0], controls[3], anc[1])
        qc.barrier()
    
    RT4L_inv(self, controls[0], controls[1], controls[2], anc[0])
    qc.barrier()
    
    if (n-4)%2 == 0:
        a_idx = 1
        for i in (range(floor((n-4)/2))):
            RT4S_inv(self, anc[a_idx - 1 + i], controls[2*i+3], controls[2*i+4], anc[a_idx + i])
            qc.barrier()
    else:
        a_idx = 2
        RTS_inv(self, anc[0], controls[3], anc[1])
        qc.barrier()
        for i in (range(floor((n-4)/2))):
            RT4S_inv(self, anc[a_idx - 1 + i], controls[2*i+4], controls[2*i+5], anc[a_idx + i])
            qc.barrier()

def apply_mct(circuit, controls, target, anc, mode='clean-ancilla'):
    if len(controls) == 1:
        circuit.cx(controls[0], target)
    elif len(controls) == 2:
        circuit.ccx(controls[0], controls[1], target)
    else:
        if mode == 'clean-ancilla':
            apply_mct_clean(circuit, controls, target, anc)
        if mode == 'dirty-ancilla':
            apply_mct_dirty(circuit, controls, target, anc)

def _mct(self, controls, target, ancilla, mode):
    if controls is None or len(controls) < 0:
        raise ValueError('you should pass controls as list or registers')
    if target is None or (not isinstance(target, Qubit) and len(target) != 1):
        raise ValueError('target length is not 1')
    if ancilla is None or len(ancilla) < ceil((len(controls)-2)/2):
        raise ValueError('for {} controls, need at least {} ancilla'.format(len(controls), ceil((len(controls)-2)/2)))
    if mode is None or (mode != 'dirty-ancilla' and mode != 'clean-ancilla'):
        raise ValueError('unknown mode')
    apply_mct(self, controls, target, ancilla, mode)

QuantumCircuit.mct = _mct
