
import os
import time
import numpy as np
import qiskit
from qiskit import QuantumCircuit
from sympy.combinatorics import Permutation
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import SabreSwap, SabreLayout, ApplyLayout, FullAncillaAllocation, EnlargeWithAncilla
import networkx as nx

import argparse
import architectures

def get_sabre_initial_map_and_swap_count(file_name, coupling_map):
    circ = QuantumCircuit.from_qasm_file(file_name)
    cm = CouplingMap(np.argwhere(coupling_map > 0))
    sabre_layout = SabreLayout(cm, seed=42)
    pm = PassManager([sabre_layout])
    mapped_circ = pm.run(circ)
    i_map =  ({circ.find_bit(q)[0] : p for (q,p) in pm.property_set['layout'].get_virtual_bits().items() if q.register.name == 'q'})
    #rolling_map  ={(p,0) : q for q,p in i_map.items()}
    rolling_map = {(p, 0) : -1 for p in range(len(coupling_map))}
    rolling_map.update({(p,0) : q for q,p in i_map.items()} )
    if 'swap' not in mapped_circ.count_ops().keys():
        swap_count = 0
    else: 
        swap_count = mapped_circ.count_ops()['swap']
    swaps = []
    cnot_counter = 0
    cnots = []
    count_dict = {0 : 0}
    for ins, qubits, clbits in mapped_circ:
        if ins.name=='cx':
            q1 = mapped_circ.find_bit(qubits[0])[0] 
            q2 = mapped_circ.find_bit(qubits[1])[0]
            cnots.append([rolling_map[q1, cnot_counter],  rolling_map[q2, cnot_counter]])
            cnot_counter += 1
            count_dict[cnot_counter] = 0
            rolling_map.update({ (p,cnot_counter) : q for (p, k),q in rolling_map.items() if k == cnot_counter-1 })
        elif ins.name=='swap':
            q1 = mapped_circ.find_bit(qubits[0])[0] 
            q2 = mapped_circ.find_bit(qubits[1])[0]
            log1 = rolling_map[q1, cnot_counter]  
            rolling_map[q1, cnot_counter]  = rolling_map[q2, cnot_counter]  
            rolling_map[q2, cnot_counter]  = log1
            
            swaps.append((q1, q2, cnot_counter, count_dict[cnot_counter]))
            count_dict[cnot_counter]  += 1
    max_per = max(count_dict.values())
    return rolling_map, swaps, cnots, max_per


def run_sabre(file_name, coupling_map):
    circ = QuantumCircuit.from_qasm_file(file_name)
    cm = CouplingMap(np.argwhere(coupling_map > 0))
    G = nx.from_numpy_array(coupling_map)

    sabre_layout = SabreLayout(cm, seed=42)
    pm = PassManager([sabre_layout])

    t_s = time.time()

    mapped_circ = pm.run(circ)
    pattern = mapped_circ.layout.routing_permutation()

    for idx0, idx1 in Permutation(pattern).transpositions():
        path = nx.dijkstra_path(G, idx0, idx1)
        for i in range(len(path)-2):
            mapped_circ.swap(path[i], path[i+1])
        mapped_circ.swap(path[-2], path[-1])
        for i in reversed(range(len(path)-2)):
            mapped_circ.swap(path[i], path[i+1])

    t_f = time.time()

    print(mapped_circ.count_ops())

    if 'swap' not in mapped_circ.count_ops().keys():
        swap_count = 0
    else: 
        swap_count = mapped_circ.count_ops()['swap']
    
    return ({'circ' : os.path.basename(file_name), 'cnots' : circ.num_nonlocal_gates(), 'swaps' : swap_count, "solve_time" : t_f - t_s, 'algorithm ': 'sabre'}, mapped_circ)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("prog", help="path to input program file")
    parser.add_argument("-o_p", "--output_path", help="where to write the resulting qasm")
    parser.add_argument("-a", "--arch", help="name of qc architecture")
    parser.add_argument('-m', '--row', type=int)
    parser.add_argument('-n', '--col', type=int)
    # parser.add_argument('-r', '--runSabre', action='store_true', help='Run Sabre')

    args = parser.parse_args()

    archs =  {
        'square': architectures.square(args.row, args.col),
        'hex': architectures.hexagon(args.row, args.col),
        'h_square': architectures.heavy_square(args.row, args.col),
        'h_hex': architectures.heavy_hexagon(args.row, args.col),
        'oct': architectures.octagon(args.row, args.col),
        'tokyo': architectures.tokyo(args.row, args.col),
    }

    if args.arch in archs:
        arch = archs[args.arch]
    base, _ = os.path.splitext(os.path.basename(args.prog))
    
    # print(arch)
    (stats, qc) = run_sabre(args.prog, arch)

    print("\n# Physical Qubit", len(arch))
    print("# Extra CNOT Count:", stats['swaps'] * 3)
    print("Total Depth:", qc.depth())

    if args.arch in archs:
        stats["arch"] = args.arch
    else:
        stats['arch'] = f"custom arch with {len(arch)} qubits"

    with open("output/" + f"stats_{base}.txt", "w") as f:
        f.write(str(stats))
        f.write("\n# Physical Qubit: " + str(len(arch)))
        f.write("\n# Extra CNOT Count: " + str(stats['swaps'] * 3))
        f.write("\nTotal Depth: " + str(qc.depth()))


    