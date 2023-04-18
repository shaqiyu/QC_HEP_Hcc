from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import ParameterVector
import numpy as np

#TODO: We need to fix the degree, it's not used for now
def FeatureMap(num_qubits=5, depth=1, degree=1, entanglement='Full', inverse=False):
     
    _num_qubits = num_qubits
    _degree = degree
    _depth = depth
    
  
    print("Checking the number of qubits:",_num_qubits)   
 
    x = ParameterVector('x',_num_qubits) 
    #x = [-9.82270658e-01, -9.65333462e-01, -9.76610780e-01, -4.06093776e-01, -4.01525736e-01, 6.00495338e-01, 5.49438477e-01]
    qr = QuantumRegister(size=_num_qubits, name='q', bits=None)
    

    qc = QuantumCircuit(qr)
    #qc = qc.bind_parameters({p: [0, 1, 2, 3, 4]})
    entangler_map = []
    for i in range(0, _num_qubits - 1, 2):
        entangler_map.append([i, i + 1])
        #print(entangler_map)
    for i in range(0, _num_qubits - 2, 2):
        entangler_map.append([i + 1, i + 2])
        #print(entangler_map)

    #if 'Full' in entanglement:
    #    for i in range(0, num_qubits, 1):
    #        for j in range(0, num_qubits, 1):
    #            if i==j:
    #                continue
    #            else:
    #                #print([i, j])
    #                if 'forward' in entanglement:
    #                    entangler_map.append([i, j])
    #                elif 'backward' in entanglement:
    #                    entangler_map.append([j, i])

    #elif 'Par' in entanglement:
    #    if 'forward' in entanglement:
    #        for i in range(0, _num_qubits, 2):
    #            entangler_map.append([i, i + 1])
    #            print(entangler_map)
    #        for i in range(0, _num_qubits - 2, 2):
    #            entangler_map.append([i + 1, i + 2]) 
    #            print(entangler_map)
    #    elif 'backward' in entanglement:
    #        for i in range(0, _num_qubits, 2):
    #            entangler_map.append([i + 1, i])
    #        for i in range(0, _num_qubits - 2, 2):
    #            entangler_map.append([i + 2, i + 1])

    #elif 'Second' in entanglement:
    #    if 'forward' in entanglement:
    #        for i in range(0, _num_qubits, 2):
    #            entangler_map.append([i, i + 1])
    #        for i in range(0, _num_qubits - 2, 2):
    #            entangler_map.append([i, i + 2])
    #    elif 'backward' in entanglement:
    #        for i in range(0, _num_qubits, 2):
    #            entangler_map.append([i + 1, i])
    #        for i in range(0, _num_qubits - 2, 2):
    #            entangler_map.append([i + 2, i])


    for d in range(_depth):
        for i in range(_num_qubits):
            #qc.t(qr[i])
            qc.h(qr[i])
            qc.rz(2*x[i], qr[i])
            #qc.u(np.pi/2.,np.pi/2., x[i], qr[i])
            #qc.t(qr[i])
            qc.ry(x[i], qr[i])
            
            #if _degree>0:
            #    qc.ry(x[i],qr[i])
        for src, targ in entangler_map:
            print(entangler_map)
            qc.cx(qr[src], qr[targ])
            #qc.sxdg(x[src], x[targ])
            #qc.rz((((x[src] + x[targ])/2)), qr[targ])
            #qc.cx(qr[src], qr[targ])
            
    if inverse:
        qc = qc.inverse()
    
    return qc

def FeatureMap_2(num_qubits=4, depth=1, degree=1, entanglement='Full', inverse=False):
    
    _num_qubits = num_qubits
    _degree = degree
    _depth = depth

    if _num_qubits > 4:
       print(" The number of qubit should be exactly for but you provided ",_num_qubits)
   
    # The vector parameters should have the same number as variables 
    x = ParameterVector('x',_num_qubits*2) 
    qr = QuantumRegister(size=_num_qubits, name='q', bits=None)

    qc = QuantumCircuit(qr)
    entangler_map = []

    if 'Full' in entanglement:
        for i in range(0, num_qubits, 1):
            for j in range(0, num_qubits, 1):
                if i==j:
                    continue
                else:    
                    #print([i, j])
                    if 'forward' in entanglement:
                        entangler_map.append([i, j])
                    elif 'backward' in entanglement:
                        entangler_map.append([j, i])

    elif 'Par' in entanglement:
        if 'forward' in entanglement:
            for i in range(0, _num_qubits, 2):
                entangler_map.append([i, i + 1])
            for i in range(0, _num_qubits - 2, 2):
                entangler_map.append([i + 1, i + 2])
        elif 'backward' in entanglement:
            for i in range(0, _num_qubits, 2):
                entangler_map.append([i + 1, i])
            for i in range(0, _num_qubits - 2, 2):
                entangler_map.append([i + 2, i + 1])

    elif 'Second' in entanglement:
        if 'forward' in entanglement:
            for i in range(0, _num_qubits, 2):
                entangler_map.append([i, i + 1])
            for i in range(0, _num_qubits - 2, 2):
                entangler_map.append([i, i + 2])
        elif 'backward' in entanglement:
            for i in range(0, _num_qubits, 2):
                entangler_map.append([i + 1, i])
            for i in range(0, _num_qubits - 2, 2):
                entangler_map.append([i + 2, i])


    for d in range(_depth):
        for i in range(_num_qubits):
            qc.h(qr[i])
            qc.rz(x[i*2], qr[i])
            #qc.t(qr[i])
            qc.ry(x[i*2+1], qr[i])
            #if _degree>0:
                #qc.ry(x[i],qr[i])
        for src, targ in entangler_map:
            qc.cx(qr[src], qr[targ])
            #qc.rx( ((x[src] + x[targ])/2)*_degree , qr[targ])
            #qc.cx(qr[src], qr[targ])
            
    if inverse:
        qc = qc.inverse()
    
    return qc
