import os
import argparse
import numpy as np
import json
from StabCode import extract, surface_code
from Stage1_Map_ver3 import S1_transpile
from Stage2_Schedule_ver3 import S2_transpile, toQASM
from qiskit import QuantumCircuit
from qiskit.visualization import plot_coupling_map
from qiskit.transpiler import CouplingMap


def coupling_to_adjacency(coupling_map):
    # Extract unique qubits and sort them
    qubits = sorted({q for pair in coupling_map for q in pair})
    qubit_index = {q: i for i, q in enumerate(qubits)}
    # Initialize adjacency matrix with zeros
    size = len(qubits)
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    # Fill in the adjacency matrix
    for a, b in coupling_map:
        i, j = qubit_index[a], qubit_index[b]
        matrix[i][j] = 1
        matrix[j][i] = 1 # Assuming undirected graph

    return matrix


def idx_dict_to_list(mapping_dict):
    # Step 1: Group by row
    from collections import defaultdict

    row_groups = defaultdict(list)
    for (row, col), val in mapping_dict.items():
        row_groups[row].append(col)

    # Step 2: Sort columns within each row
    for row in row_groups:
        row_groups[row].sort()

    # Step 3: Convert to [row_index, sequential_index] list
    result = []
    for row in sorted(row_groups.keys()):
        for seq_index, col in enumerate(row_groups[row]):
            result.append([row - min(row_groups.keys()), seq_index])

    return result


def transpile_circuit_QECCSynth(code_distance, architecture):
    """Generate a qiskit Circuit given a surface code (distance) and a hardware architecture (coupling graph)

    Args:
        code_distance (int): Distance of the surface code
        backend (Tuple): Tuple containing coupling graph and index2coordinates dictionary

    Returns:
        QuantumCircuit: Qiskit QuantumCircuit object with constructed circuit
    """

    protocol = 'Flag-Bridge'
    
    maxL = None # maximum size for the ancilla block (If none, then this value is computed)
    maxT = 50 # upper bound for time step (minutes or seconds?)
    Len = None # not defined what this is ...
    chunkNum = 1 # Not defined what this is...


    # Generate surface code (written as stabilizers to a file), given code distance
    
    #code_distance = 3
    stabilizer_file = f'surface_code_distance_{code_distance}'
    codeName = f'surface_{code_distance}'
    surface_code(stabilizer_file, code_distance)
    stabilizers, dataNum = extract(stabilizer_file)
    code = (dataNum, stabilizers)
    #print(code)
    


    #stabilizers, dataNum = extract("../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/code/Surface_3.stab")
    #code = (dataNum, stabilizers)
    #print(code)


    
    # TODO: generate arch from backend (Coupling graph)
    archName = "square_chiplets"

    # Coupling graph and Qubit mapping
    CG, idx2coord = architecture

    # Plot
    #cm = CouplingMap(couplinglist=CG)
    #plot_coupling_map(cm.size(), None, cm.get_edges(), filename="testing.png")
    
    CG = np.array(coupling_to_adjacency(CG))
    idx2coord = idx_dict_to_list(idx2coord)

    #print(CG)
    #print(idx2coord_mech)

    #with open("../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/arch/square_17x17.json", 'r') as f:
    #    arch = json.load(f)
    #    CG = np.array(arch['CG'])
    #    idx2coord = arch['idx2coord']

    #print(idx2coord)
    #print(CG)
    #print(idx2coord)
    #return None
    base = codeName + '-' + archName

    physNum = len(CG)
    print("\n# Physical Qubit: %d" % physNum)

    # 1. Stage
    print('\n\nStage 1:')
    S1_transpile(code, CG, protocol, maxL, Len, chunkNum, '../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/aux_files/clause_'+base, '../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/aux_files/sol_'+base, '../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/output/S1_'+base)
    with open('../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/output/S1_'+base+'.json', 'r') as f:
        result = json.load(f)

    physNum = len(CG)
    print("\n# Physical Qubit: %d" % physNum)

    # 2. Stage
    print('\n\nStage 2:')
    S2_transpile(result['chunkNum'], result['Data'], result['Stab'], CG, idx2coord, result['Collison_list'], maxT, '../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/output/S2_'+base)

    with open('../qecc_mapping/src/baseline/QECC_Synth/SurfStitch/MyCode/output/S2_'+base+'.json', 'r') as f:
        result = json.load(f)    

    print('\nResult:')
    print('\nData Mapping:')
    for data in range(dataNum):
        print('q[%s]' % data, end='\t')
        for k in range(result['chunkNum']):
            print('-> %d' % (result['Data'][k][data]), end='\t')
        print()

    cnotNum = 0

    print('\nSwap Layer:')
    for k in range(result['chunkNum']):
        # print('Chunk %d <-> Chunk %d:' % (k, (k+1)%result['chunkNum']))
        for s in result['Swap_layer'][k]:
            # print('SWAP(%d, %d)' % (s[1][1], s[1][2]), end='\t')
            cnotNum += 3
        print()

    print('\nBridge Allocating: Stabilizer -> Syndrome qubit -> Ancilla bridge')
    for k in range(result['chunkNum']):
        for k, stab in enumerate(result['Stab'][k]):
            print('%s -> %s -> %s' % (stab['name'], stab['synd'], set(stab['Ancilla'])))
            cnotNum += len(stab['Ancilla']) * 2 - 2

    physNum = len(CG)
    print("\n# Physical Qubit: %d" % physNum)
    print("Code Density: %.2f" % (2*(sum(len(stab['Ctrl']) for stab in stabilizers))/(dataNum+len(stabilizers))))
    print('# Extra Cnot Gates: %d' % cnotNum)

    QC = QuantumCircuit(physNum, physNum)
    for i in range(result['chunkNum']):
        QC = QC.compose(toQASM(len(CG), result['Circuit'][i]))
        QC = QC.compose(toQASM(len(CG), result['Swap_layer'][i]))

    print('Circuit Depth: %d' % QC.depth())
    print()


    # Save transpiled circuit
    import pickle

    filehandler = open('transpiled_circuit_qecc_synth', 'wb') 
    pickle.dump(QC, filehandler)
    print("pickled file!")


    return QC