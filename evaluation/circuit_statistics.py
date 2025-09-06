import sys
import os
import qiskit
import networkx as nx

# MECH
sys.path.append(os.path.join(os.getcwd(), "./src/baseline/MECH"))
from src.baseline.MECH import Router
from src.baseline.MECH.Chiplet import *
from src.baseline.MECH.Circuit import *


def calc_circuit_qiskit_stats(transpiled_circuit: qiskit.circuit, G: nx.graph) -> dict:
    """_summary_

    Args:
        transpiled_circuit (qiskit.circuit): _description_
        G (nx.graph): _description_

    Returns:
        dict: _description_
    """

    # Parameter
    cross_chip_gate_weight = 7.4

    # Decompose swap gates
    filter_function = lambda gate: gate.operation.num_qubits >= 2
    swap_decomposed_circuit = transpiled_circuit.decompose('swap')

    # Calculate depth
    swap_decomposed_depth = 0
    swap_decomposed_depth += swap_decomposed_circuit.depth(filter_function)

    # Count instructions
    for instr, qargs, cargs in swap_decomposed_circuit.data:
        # Only care about 2-qubit gates
        if instr.num_qubits < 2:
            continue

        # Extract physical qubit indices
        q1 = swap_decomposed_circuit.qubits.index(qargs[0])
        q2 = swap_decomposed_circuit.qubits.index(qargs[1])

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

    # Calculate effective number of CNOT gates (for calculation, see section 7.1 of "MECH: Multi-Entry Communication Highway for Superconducting Quantum Chiplets")
    norm_cnots = within_chip_cnots + cross_chip_cnots * cross_chip_gate_weight

    # Collect statistics
    # print('Qiskit: decomposed_depth = {}, within_chip_cnots={}, cross_chip_cnots={}, norm_cnots = {}'.format(swap_decomposed_depth, within_chip_cnots, cross_chip_cnots, norm_cnots))
    result_qiskit = {'depth': swap_decomposed_depth, 'eff_gate_num': norm_cnots, 'on-chip': within_chip_cnots, 'cross-chip': cross_chip_cnots}
    
    return result_qiskit


def calc_circuit_mech_stats(router: Router) -> dict:
    """_summary_

    Args:
        router (Router): _description_

    Returns:
        dict: _description_
    """
    
    # Parameter
    cross_chip_gate_weight = 7.4
    meas_weight = 2.2

    on_chip_gate_num = 0
    cross_chip_gate_num = 0
    meas_num = 0

    for idx in range(router.circuit.depth):
        for line in range(len(router.circuit.circuit_lines)):
            if router.circuit.take_role(line, idx) == 'q':
                meas_num += 1
            if router.circuit.take_role(line, idx) in ['t', 'mt']:
                node = router.circuit.take_node(line, idx)
                if isinstance(node, OpNode):
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

    # Collect statistics
    # print('MECH: decomposed_depth = {}, within_chip_cnots={}, cross_chip_cnots={}, norm_cnots = {}'.format(router.circuit.depth, on_chip_gate_num, cross_chip_gate_num, eff_gate_num))
    result_mech = {'depth': router.circuit.depth, 'eff_gate_num': eff_gate_num, 'on-chip': on_chip_gate_num, 'cross-chip': cross_chip_gate_num, 'meas_num': meas_num, 'shuttle_num': len(router.highway_manager.shuttle_stack)}

    return result_mech