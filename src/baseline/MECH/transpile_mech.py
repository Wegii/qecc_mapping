import qiskit.circuit
from Circuit import *
from Chiplet import *
from HighwayOccupancy import *
from Router import *
from MECHBenchmarks import *
import json
import os

from tqdm import tqdm
import pickle

import matplotlib.pyplot as plt

import qiskit



def QFT(data_qubit_num):
    program = []
    for i in range(data_qubit_num-1):
        control_block = ControlBlock(i, Counter(list(range(i+1, data_qubit_num))*2))
        program.append(control_block)
    
    return program


def generate_simple_chiplet():
    """
    Generate simplified backend with highway layout and chiplet connections
    """
    
    lane_num = 1
    structure = 'square'
    chip_col_num = 2
    chip_row_num = 2
    x_num, y_num = 9,9
    # Generate layout of chip, by connecting smaller chiplets
    G = gen_chiplet_array(structure, chip_col_num, chip_row_num, x_num, y_num, cross_link_sparsity=1)
    # Add highway
    gen_highway_layout(G)
    return G


def circuit_to_program(circuit):
    """ Circuit to program 
    
    TODO: remove data_qubit_num, as this is generally not necessary
    """

    program = []
    for instr, qargs, cargs in circuit.data:
        
        if len(qargs) == 2:
            # Add two-qubit gate    
            control = circuit.qubits.index(qargs[0])
            target = circuit.qubits.index(qargs[1])
            program.append(ControlBlock(control, Counter(list([target]))))
    #print(len(program))

    return program


def program_to_circuit(router):
    """_summary_

    Args:
        router (_type_): _description_

    Returns:
        _type_: _description_
    """

    # Construct circuit with correct number of qubits
    circuit = qiskit.QuantumCircuit(len(router.circuit.circuit_lines))

    # Iterate over depth
    for idx in range(router.circuit.depth):
        # Iterate over each qubit
        for line in range(len(router.circuit.circuit_lines)):
            #print(router.circuit.take_role(line, idx))
            
            if router.circuit.take_role(line, idx) == 'q':
                #meas_num += 1
                # Single qubit-gate
                node = router.circuit.take_node(line, idx)
                control_line = node.q
                control_qubit = router.highway_manager.idx_qubit_dict[control_line]
                # Measurement in X-basis
                circuit.H(control_qubit)
                circuit.reset(control_qubit)

            # Iterate over targets
            """
            if router.circuit.take_role(line, idx) in ['t', 'mt']:
                # Target qubit

                node = router.circuit.take_node(line, idx)
                if isinstance(node,  OpNode):
                    control_line = node.control
                if isinstance(node, MOpNode):
                    control_line = node.shared
                control_qubit, target_qubit = router.highway_manager.idx_qubit_dict[control_line], router.highway_manager.idx_qubit_dict[line]

                circuit.ControlledGate(control_qubit, target_qubit)
            """
    
    #print(on_chip_gate_num)
    print("Conversion succeeded")
    return circuit


def transpile_circuit_MECH(circuit, backend=None):
    """ Transpile using MECH approach

    backend = backend to use; e. g. nighthawk. Currently not used
    """
    
    transpiled_qc = circuit


    # Generate backend with highway
    G = generate_simple_chiplet()
    gen_highway_layout(G)

    # Translate circuit to program (necessary for MECH)
    #program = QFT(len(G.nodes) - len(get_highway_qubits(G))) # circuit_to_program(circuit)
    program = circuit_to_program(circuit)
    print("Program generated")
    
    
    prep_period = 13
    meas_period = 2
    cross_chip_gate_weight = 7.4
    router = Router(G, prep_period=prep_period, meas_period=meas_period)


    if False:
        for i, control_block in enumerate(tqdm(program, total=len(program))):

            #if target != control
            router.execute_control_multi_target_block(control_block, cross_chip_gate_weight, cross_chip_overhead = 3)
            
            #else:
            """
                earliest_idx = self.get_highway_aware_earliest_index_for_1qubit_op(circuit, op)
                if earliest_idx < 0:
                    earliest_idx = self.get_shuttle_exec_start_time(len(self.shuttle_stack)-1)
                circuit.add_1qubit_op(op, depth=earliest_idx + 1)

                op = OpNode('M')
                router.circuit.add_1qubit_op(op, depth=None)
                
                def add_1qubit_op(self, op, depth=None):
                    if depth is not None:
                        idx = depth - 1
                    else:
                        idx = self.get_line_depth(op.q)
                    self.add_node_with_role(op.q, idx, op, 'q')
                
            """
        #print(router.circuit.circuit_lines)
        #print(len(router.highway_manager.shuttle_stack))


        for shuttle_idx in range(len(router.highway_manager.shuttle_stack)):
            for idx in range(router.highway_manager.get_shuttle_prep_start_time(shuttle_idx), router.highway_manager.get_shuttle_exec_start_time(shuttle_idx)):
                for line in range(len(router.circuit.circuit_lines)):
                    assert router.circuit.is_position_empty(line, idx)

            source_measured_qubits_dict = router.highway_manager.bridge_throughout_highway(router.circuit, shuttle_idx)
            router.highway_manager.reentangle_throughout_highway(router.circuit, shuttle_idx, source_measured_qubits_dict)
        

        for shuttle_idx in range(len(router.highway_manager.shuttle_stack)):
            router.highway_manager.measure_throughout_highway(router.circuit, shuttle_idx)



    with open('../qecc_mapping/tests/transpiled_circuit_mech', 'rb') as file:
        router_loaded = pickle.load(file)

    print("loaded router")
    router = router_loaded
    #transpiled_qc = program_to_circuit(router)
    transpiled_qc_testing = program_to_circuit(router)
    #print(transpiled_qc_testing)


    #filehandler = open('transpiled_circuit_mech', 'wb') 
    #pickle.dump(router, filehandler)
    #print("pickled file")

    return transpiled_qc