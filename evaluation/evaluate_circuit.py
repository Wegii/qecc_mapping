import sys
import os
sys.path.append(os.path.join(os.getcwd(), "."))

import pickle
import itertools
from tqdm import tqdm
import logging

# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
from circuit_statistics import *
from experiments.transpilation_experiments import generate_simple_dqc_backend

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


def plot_gate_overhead(result_mech: dict, result_qecc_synth: dict, result_qiskit: dict, filename) -> None:

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
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


base_path = "./data/transpiled_circuit/"
code_name = "surface"

# Iterate over different backend configurations
nnx = [4, 6, 8]
cl = [1, 2]
for n, c in tqdm(itertools.product(nnx, cl), total=len(nnx) * len(cl)):
    # TODO: load files based on backend configuration

    # Load architecture
    simple_dqc_backend, qubit_num, data_qubit_num = generate_simple_dqc_backend(n, n, c)
    
    # MECH
    logging.info("MECH: Loading file")
    with open(base_path + f'{code_name}_d_square_{n}_{n}_{c}_mech', 'rb') as file:
        circuit_mech = pickle.load(file)
    logging.info("MECH: Calculating statistics")
    result_mech = calc_circuit_mech_stats(circuit_mech)

    #QECC-Synth
    #with open(base_path + f'{code_name}_d_square_{n}_{n}_{c}_qecc_synth', 'rb') as file:
    #    circuit_mech = pickle.load(file)

    # SABRE
    logging.info("SABRE: Loading file")
    with open(base_path + f'{code_name}_d_square_{n}_{n}_{c}_qiskit', 'rb') as file:
        circuit_qiskit = pickle.load(file)
    logging.info("SABRE: Calculating statistics")
    result_qiskit = calc_circuit_qiskit_stats(circuit_qiskit, simple_dqc_backend)
    
    # Plot results
    filename = f'evaluation/figures/gate_overhead_{code_name}_d_square_{n}_{n}_{c}.png'
    plot_gate_overhead(result_mech, None, result_qiskit, filename)
    