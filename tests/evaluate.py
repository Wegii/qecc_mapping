# Hacky way of adding project. Change this at some point
import os
import sys
sys.path.append(os.path.join(os.getcwd(), "qecc_mapping"))
sys.path.append(os.path.join(os.getcwd(), "qecc_mapping/src/baseline/MECH"))

# Include MECH transpilation
from src.baseline.MECH.Circuit import *
from src.baseline.MECH.Chiplet import *
from src.baseline.MECH.HighwayOccupancy import *
from src.baseline.MECH.Router import *
from src.baseline.MECH.MECHBenchmarks import *

import pickle

# eccentric_bench
sys.path.append(os.path.join(os.getcwd(), "eccentric_bench"))
sys.path.append(os.path.join(os.getcwd(), "eccentric_bench/utils"))
sys.path.append(os.path.join(os.getcwd(), "eccentric_bench/utils/math.py"))
sys.path.append(os.path.join(os.getcwd(), "eccentric_bench/external/qiskit_qec/src"))
from utils.math import *
from qiskit.compiler import transpile
from codes import get_code, get_max_d
from transpilers import translate

# Qiskit
from qiskit.transpiler import CouplingMap


# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns


marker_styles = {
    'gross': 'o',
    'bacon': 's',
    'hh': '^',
    'surface': 'D',
    'other': 'X'
}

code_rename_map = {
    'bacon': 'Bacon-Shor',
    'hh': 'Heavy-hex',
    'gross': 'Gross'
}


def generate_qiskit_surface_code_circuit(mech_backend):
    code_name = "surface"
    translating_method = None
    
    # Generate qec code
    data_qubit_num = len(mech_backend.nodes) - len(get_highway_qubits(mech_backend))
    d = get_max_d(code_name, data_qubit_num)
    cycles = d
    code = get_code(code_name, d, cycles)

    if translating_method:
        code.qc = translate(code.qc, translating_method)
        #TODO: either else here or sth

    original_circuit = code.qc

    return original_circuit


def generate_qiskit_backend_from_mech():
    structure = 'square'
    chip_col_num = 2
    chip_row_num = 2
    x_num, y_num = 9,9
    # Generate layout of chip, by connecting smaller chiplets
    G = gen_chiplet_array(structure, chip_col_num, chip_row_num, x_num, y_num, cross_link_sparsity=1)
    # Add highway
    gen_highway_layout(G)

    qubit_idx_dict = gen_qubit_idx_dict(G)
    regular_coupling = list([qubit_idx_dict[n1], qubit_idx_dict[n2]] for n1,n2 in G.edges)
    regular_coupling += list([qubit_idx_dict[n2], qubit_idx_dict[n1]] for n1,n2 in G.edges)

    return regular_coupling, G


def calc_circuit_stats(router, transpiled_circuit, mech_backend):
    
    cross_chip_gate_weight = 7.4
    meas_weight = 2.2

    # MECH
    on_chip_gate_num = 0
    cross_chip_gate_num = 0
    meas_num = 0

    for idx in range(router.circuit.depth):
        for line in range(len(router.circuit.circuit_lines)):
            if router.circuit.take_role(line, idx) == 'q':
                meas_num += 1
            if router.circuit.take_role(line, idx) in ['t', 'mt']:
                node = router.circuit.take_node(line, idx)
                if isinstance(node,  OpNode):
                    control_line = node.control
                if isinstance(node, MOpNode):
                    control_line = node.shared
                control_qubit, target_qubit = router.highway_manager.idx_qubit_dict[control_line], router.highway_manager.idx_qubit_dict[line]
                
                if router.chip.has_edge(control_qubit, target_qubit) and router.chip.edges[(control_qubit, target_qubit)]['type'] == 'cross_chip':
                    cross_chip_gate_num += 1
                else:
                    on_chip_gate_num += 1

    # Calculate effective number of CNOT gates (for calculation, see section 7.1 of "MECH: Multi-Entry Communication Highway for Superconducting Quantum Chiplets")
    eff_gate_num = on_chip_gate_num + cross_chip_gate_num * cross_chip_gate_weight + meas_num * meas_weight
    print(f"Eff_CNOT: {eff_gate_num}")

    # Collect statistics
    result_mech = {'depth': router.circuit.depth, 'eff_gate_num': eff_gate_num, 'on-chip': on_chip_gate_num, 'cross-chip': cross_chip_gate_num, 'meas_num': meas_num, 'shuttle_num': len(router.highway_manager.shuttle_stack)}
    print('MECH: decomposed_depth = {}, within_chip_cnots={}, cross_chip_cnots={}, norm_cnots = {}'.format(router.circuit.depth, on_chip_gate_num, cross_chip_gate_num, eff_gate_num))

    # Qiskit
    transpiled_depth = 0
    swap_decomposed_depth = 0
    within_chip_cnots = 0
    cross_chip_cnots = 0
    norm_cnots = 0
    filter_function = lambda gate: gate.operation.num_qubits >= 2

    transpiled_depth += transpiled_circuit.depth()
    swap_decomposed_circuit = transpiled_circuit.decompose('swap')
    swap_decomposed_depth += swap_decomposed_circuit.depth(filter_function)

    def count_norm_cnots(G, sabre_circ, cross_chip_ratio=7.4):
        within_chip_cnots = 0
        cross_chip_cnots = 0
        for instr, qargs, cargs in sabre_circ.data:
            # Only care about 2-qubit gates
            if instr.num_qubits < 2:
                continue

            #print(instr, qargs)

            # Extract physical qubit indices
            #q1, q2 = qargs[0].index, qargs[1].index
            q1 = sabre_circ.qubits.index(qargs[0])
            q2 = sabre_circ.qubits.index(qargs[1])

            idx_qubit_dict = gen_idx_qubit_dict(G)
            edge = (idx_qubit_dict[q1], idx_qubit_dict[q2])

            if instr.name in ['cx', 'cp'] and G.edges[edge]['type'] == 'on_chip':
                within_chip_cnots += 1
            elif instr.name == 'swap' and G.edges[edge]['type'] == 'on_chip':
                within_chip_cnots += 3
            elif instr.name in ['cx', 'cp'] and G.edges[edge]['type'] == 'cross_chip':
                cross_chip_cnots += 1
            elif instr.name == 'swap' and G.edges[edge]['type'] == 'cross_chip':
                cross_chip_cnots += 3
                
        return within_chip_cnots, cross_chip_cnots, within_chip_cnots + cross_chip_cnots * cross_chip_ratio
    
    w, c, n = count_norm_cnots(mech_backend, swap_decomposed_circuit, cross_chip_gate_weight)
    within_chip_cnots += w
    cross_chip_cnots += c
    norm_cnots += n

    result_qiskit = {'depth': swap_decomposed_depth, 'eff_gate_num': norm_cnots, 'on-chip': within_chip_cnots, 'cross-chip': cross_chip_cnots}
    print('Qiskit: decomposed_depth = {}, within_chip_cnots={}, cross_chip_cnots={}, norm_cnots = {}'.format(swap_decomposed_depth, within_chip_cnots, cross_chip_cnots, norm_cnots))

    return result_mech, result_qiskit


def plot_gate_overhead(result_mech=None, result_qiskit=None):
    # Load data
    codes = ['MECH', 'Qiskit-SABRE']
    # On-chip CNOT, cross-chip CNOT, eff_cnot, depth
    mech_vals = [result_mech['on-chip'], result_mech['cross-chip'], result_mech['eff_gate_num'], result_mech['depth']]
    sabre_vals = [result_qiskit['on-chip'], result_qiskit['cross-chip'], result_qiskit['eff_gate_num'], result_qiskit['depth']]

    metrics = ["On-chip CNOT", "Cross-chip CNOT", "Eff_CNOT", "Depth"]

    x = np.arange(len(codes))  # two groups (MECH, SABRE)
    n_metrics = len(metrics)
    width = 0.8 / n_metrics   # distribute 4 bars in each group

    # Colors & hatches
    base_palette = sns.color_palette("pastel", n_colors=n_metrics)
    hatches = ["/", "\\", "o", "x"]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle("Code: Surface", fontsize=16)

    # Plot bars for each metric
    for i, (metric, color, hatch) in enumerate(zip(metrics, base_palette, hatches)):
        ax.bar(
            x + (i - n_metrics/2) * width + width/2,  # offset to center bars in each group
            [mech_vals[i], sabre_vals[i]],
            width,
            label=metric,
            color=color,
            edgecolor="black",
            hatch=hatch
        )

    # X-axis
    ax.set_xticks(x)
    ax.set_xticklabels(codes)

    # Y-axis
    ax.set_ylabel("Number of gates")
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.7)

    ax.text(-0.1, 1.05, 'Lower is better ↓', transform=ax.transAxes, fontsize=10, fontweight='bold', va='top', ha='left')
    # Legend
    ax.legend(title="Metric", fontsize=12, title_fontsize=12, frameon=True, ncol=2)

    plt.tight_layout()
    plt.savefig(f"gate_overhead_mech_vs_sabre.png", bbox_inches="tight")
    plt.close()


# Qiskit
# Generate backend
regular_coupling, backend = generate_qiskit_backend_from_mech()
# Load original circuit
original_circuit = generate_qiskit_surface_code_circuit(backend)
# Transpile
transpiled_circuit = transpile(original_circuit, coupling_map=CouplingMap(regular_coupling), optimization_level=2)


# MECH
# Load transpiled circuit
with open('qecc_mapping/tests/transpiled_circuit_mech', 'rb') as file:
    router_loaded = pickle.load(file)


# Calculate statistics
result_mech, result_qiskit = calc_circuit_stats(router_loaded, transpiled_circuit, backend)
print(result_mech)

plot_gate_overhead(result_mech, result_qiskit)
#plot_gate_overhead()

print("Done!")