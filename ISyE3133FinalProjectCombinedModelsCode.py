# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 21:47:40 2024

@author: sydbi
"""

# Testing is at the end of the code
from gurobipy import Model, GRB
import gurobipy as gp
import numpy as np
import matplotlib.pyplot as plt

# Model 1 Function with input String S representing the nucelotide sequence
def solve_nucleic_acid_folding_model_1(S):
    n = len(S)

    m = Model("nucleic_acid_folding")

    # Variables: x[i, j] is 1 if nucleotide i pairs with nucleotide j
    x = {}
    for i in range(n):
        for j in range(i + 1, n):
            if (S[i] == 'A' and S[j] == 'U') or (S[i] == 'U' and S[j] == 'A') or \
            (S[i] == 'G' and S[j] == 'C') or (S[i] == 'C' and S[j] == 'G'):
                x[i, j] = m.addVar(vtype=GRB.BINARY, name=f"x[{i},{j}]")

    # Objective: Maximize number of stable pairings
    m.setObjective(sum(x.values()), GRB.MAXIMIZE)

    # Non-Crossing Constraint for a Circular Configuration
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in x:
                for k in range(i + 1, j):
                    for l in range(j + 1, n + i):
                        actual_l = l % n
                        if (k, actual_l) in x and actual_l > k:
                            m.addConstr(x[i, j] + x[k, actual_l] <= 1, f"non_crossing[{i},{j},{k},{actual_l}]")

    # Single Pairing Constraint
    for i in range(n):
        m.addConstr(sum(x.get((min(i, j), max(i, j)), 0) for j in range(n) if j != i) <= 1, f"single_pairing[{i}]")

    m.optimize()

    pairs = []
    num_AU_pairs = 0
    num_GC_pairs = 0

    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in x and x[i, j].X > 0.5:
                pairs.append((i, j))
                if (S[i], S[j]) == ('A', 'U') or (S[i], S[j]) == ('U', 'A'):
                    num_AU_pairs += 1
                elif (S[i], S[j]) == ('G', 'C') or (S[i], S[j]) == ('C', 'G'):
                    num_GC_pairs += 1

    # Results
    print("Model 1")
    print("the objective")
    print(m.objVal)
    print("Number of AU pairs")
    print(num_AU_pairs)
    print("Number of GC pairs")
    print(num_GC_pairs)

    # Visualization
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_plot = 3 * np.cos(theta)
    y_plot = 3 * np.sin(theta)

    x_lab = 4 * np.cos(theta)
    y_lab = 4 * np.sin(theta)

    x_ind = 3.5 * np.cos(theta)
    y_ind = 3.5 * np.sin(theta)

    fig, ax = plt.subplots(figsize=(8, 8))
    circle = plt.Circle((0, 0), 3, fill=False, linewidth=1, linestyle='--')
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()
    ax.add_patch(circle)
    ax.scatter(x_plot, y_plot, color='black')

    for i in range(n):
        ax.text(x_lab[i], y_lab[i], S[i], fontsize=12, ha='center', va='center')
        ax.text(x_ind[i], y_ind[i], str(i), ha='center', va='center', fontsize=8, color='red')

    for i, j in pairs:
        ax.plot([x_plot[i], x_plot[j]], [y_plot[i], y_plot[j]], 'g-')

    plt.show()


# Model 2 Function with input String S representing the nucelotide sequence
def solve_nucleic_acid_folding_model_2(S):
    n = len(S)
    m = Model("nucleic_acid_folding_model_2")

    # Pair weights
    weights = {
        ('G', 'C'): 3.0,
        ('A', 'U'): 2.0,
        ('G', 'U'): 0.1,
        ('A', 'C'): 0.05
    }

    # Variables
    x = {}
    for i in range(n):
        for j in range(i + 3, n):  # Enforcing distance of at least 3
            if (S[i], S[j]) in weights or (S[j], S[i]) in weights:
                x[i, j] = m.addVar(vtype=GRB.BINARY, name=f"x[{i},{j}]")

    # Objective: Maximize the weighted sum of selected pairs
    m.setObjective(sum(weights.get((S[i], S[j]), weights.get((S[j], S[i]), 0)) * x[i, j] for i, j in x), GRB.MAXIMIZE)

    # Non-Crossing Constraint for a Circular Configuration
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                for k in range(i + 1, j):
                    for l in range(j + 1, n + k):  # Including wrap-around for circularity
                        if l < n and (k, l) in x:
                            m.addConstr(x[i, j] + x[k, l] <= 1)

    # Single Pairing Constraint
    for i in range(n):
        m.addConstr(sum(x.get((min(i, j), max(i, j)), 0) for j in range(n) if i != j) <= 1)

    m.optimize()
    pairs = [(i, j) for i, j in x if x[i, j].X > 0.5]
    num_AU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'U'), ('U', 'A')])
    num_GC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'C'), ('C', 'G')])
    num_GU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'U'), ('U', 'G')])
    num_AC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'C'), ('C', 'A')])

    print("Model 2")
    print("the objective")
    print(m.objVal)
    print("Number of AU pairs")
    print(num_AU_pairs)
    print("Number of GC pairs")
    print(num_GC_pairs)
    print("Number of GU pairs")
    print(num_GU_pairs)
    print("Number of AC pairs")
    print(num_AC_pairs)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_plot = 3 * np.cos(theta)
    y_plot = 3 * np.sin(theta)

    x_lab = 4 * np.cos(theta)
    y_lab = 4 * np.sin(theta)

    x_ind = 3.5 * np.cos(theta)
    y_ind = 3.5 * np.sin(theta)

    fig, ax = plt.subplots(figsize=(8, 8))
    circle = plt.Circle((0, 0), 3, fill=False, linewidth=1, linestyle='--')
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()
    ax.add_patch(circle)
    ax.scatter(x_plot, y_plot, color='black')

    for i in range(n):
        ax.text(x_lab[i], y_lab[i], S[i], fontsize=12, ha='center', va='center')
        ax.text(x_ind[i], y_ind[i], str(i), ha='center', va='center', fontsize=8, color='red')

    for i, j in pairs:
        ax.plot([x_plot[i], x_plot[j]], [y_plot[i], y_plot[j]], 'g-')

    plt.show()
    return pairs



# Model 3 Function with input String S representing the nucelotide sequence
def solve_nucleic_acid_folding_model_3(S):
    n = len(S)
    m = Model("nucleic_acid_folding_model_3")

    # Pair weights
    weights = {
        ('G', 'C'): 3.0,
        ('A', 'U'): 2.0,
        ('G', 'U'): 0.1,
        ('A', 'C'): 0.05
    }

    # Variables
    x = {}
    for i in range(n):
        for j in range(i + 3, n):  # Enforcing distance of at least 3
            if (S[i], S[j]) in weights or (S[j], S[i]) in weights:
                x[i, j] = m.addVar(vtype=GRB.BINARY, name=f"x[{i},{j}]")

    y = {}
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                if (i + 1, j - 1) in x:
                    y[i, j] = m.addVar(vtype=GRB.BINARY, name=f"y[{i},{j}]")


    y_total = sum(y[i, j] for i in range(n) for j in range(i + 3, n) if (i, j) in y)

    # Objective: Maximize the weighted sum of selected pairs
    m.setObjective(sum(weights.get((S[i], S[j]), weights.get((S[j], S[i]), 0)) * x[i, j] for i, j in x) + y_total, GRB.MAXIMIZE)

    # Non-Crossing Constraint for a Circular Configuration
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                for k in range(i + 1, j):
                    for l in range(j + 1, n + k):  # Including wrap-around for circularity
                        if l < n and (k, l) in x:
                            m.addConstr(x[i, j] + x[k, l] <= 1)

    # Single Pairing Constraint
    for i in range(n):
        m.addConstr(sum(x.get((min(i, j), max(i, j)), 0) for j in range(n) if i != j) <= 1)

    # Stacked Quartet Constraints
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                if (i + 1, j - 1) in x:
                    m.addConstr(x[i, j] + x[i + 1, j - 1] - y[i, j] <= 1)
                    m.addConstr((2 * y[i, j]) - x[i, j] - x[i + 1, j - 1] <= 0)


    m.optimize()

    pairs = [(i, j) for i, j in x if x[i, j].X > 0.5]
    num_AU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'U'), ('U', 'A')])
    num_GC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'C'), ('C', 'G')])
    num_GU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'U'), ('U', 'G')])
    num_AC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'C'), ('C', 'A')])

    num_stacked_quartets = 0
    for (i, j) in pairs:
        if (i + 1, j - 1) in pairs:
            num_stacked_quartets += 1

    print("Model 3")
    print("the objective")
    print(m.objVal)
    print("Number of AU pairs")
    print(num_AU_pairs)
    print("Number of GC pairs")
    print(num_GC_pairs)
    print("Number of GU pairs")
    print(num_GU_pairs)
    print("Number of AC pairs")
    print(num_AC_pairs)
    print("Number of Stacked Quartets")
    print(num_stacked_quartets)


    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_plot = 3 * np.cos(theta)
    y_plot = 3 * np.sin(theta)

    x_lab = 4 * np.cos(theta)
    y_lab = 4 * np.sin(theta)

    x_ind = 3.5 * np.cos(theta)
    y_ind = 3.5 * np.sin(theta)

    fig, ax = plt.subplots(figsize=(8, 8))
    circle = plt.Circle((0, 0), 3, fill=False, linewidth=1, linestyle='--')
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()
    ax.add_patch(circle)
    ax.scatter(x_plot, y_plot, color='black')

    for i in range(n):
        ax.text(x_lab[i], y_lab[i], S[i], fontsize=12, ha='center', va='center')
        ax.text(x_ind[i], y_ind[i], str(i), ha='center', va='center', fontsize=8, color='red')

    for i, j in pairs:
        ax.plot([x_plot[i], x_plot[j]], [y_plot[i], y_plot[j]], 'g-')

    plt.show()
    return pairs


# Model 4 Function with input String S representing the nucelotide sequence
def solve_nucleic_acid_folding_model_4(S):
    n = len(S)
    m = Model("nucleic_acid_folding_model_4")

    # Pair weights
    weights = {
        ('G', 'C'): 3.0,
        ('A', 'U'): 2.0,
        ('G', 'U'): 0.1,
        ('A', 'C'): 0.05
    }

    # Variables
    x = {}
    for i in range(n):
        for j in range(i + 3, n):  # Enforcing distance of at least 3
            if (S[i], S[j]) in weights or (S[j], S[i]) in weights:
                x[i, j] = m.addVar(vtype=GRB.BINARY, name=f"x[{i},{j}]")

    y = {}
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                if (i + 1, j - 1) in x:
                    y[i, j] = m.addVar(vtype=GRB.BINARY, name=f"y[{i},{j}]")

    y_total = sum(y[i, j] for i in range(n) for j in range(i + 3, n) if (i, j) in y)

    c = {}
    for i in range(n):
        for j in range(n):
            if (i, j) in x:
                for k in range(i + 1, j):
                    for l in range(j + 1, n + k):
                        if l < n and (k, l) in x:
                            c[i, j, k, l] = m.addVar(vtype=GRB.BINARY, name=f"c[{i},{j},{k},{l}]")

    c_total = sum(c[i, j, k, l] for i in range(n) for j in range(n) for k in range(n) for l in range(n + k) if (i, j, k, l) in c)


    # Objective: Maximize the weighted sum of selected pairs
    m.setObjective(sum(weights.get((S[i], S[j]), weights.get((S[j], S[i]), 0)) * x[i, j] for i, j in x) + y_total, GRB.MAXIMIZE)


    # c(i, j, k, l) exists
    for i in range(n):
        for j in range(n):
            if (i, j) in x:
                for k in range(i + 1, j):
                    for l in range(j + 1, n + k):  # Including wrap-around for circularity
                        if l < n and (k, l) in x:
                            m.addConstr(x[i, j] + x[k, l] - c[i, j, k, l] <= 1)

    # Limit of 10 Crossing Pairs Constraint
    m.addConstr(c_total <= 10)

    # Single Pairing Constraint
    for i in range(n):
        m.addConstr(sum(x.get((min(i, j), max(i, j)), 0) for j in range(n) if i != j) <= 1)

    # Stacked Quartet Constraints
    for i in range(n):
        for j in range(i + 3, n):
            if (i, j) in x:
                if (i + 1, j - 1) in x:
                    m.addConstr(x[i, j] + x[i + 1, j - 1] - y[i, j] <= 1)
                    m.addConstr((2 * y[i, j]) - x[i, j] - x[i + 1, j - 1] <= 0)

    m.optimize()

    pairs = [(i, j) for i, j in x if x[i, j].X > 0.5]
    num_AU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'U'), ('U', 'A')])
    num_GC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'C'), ('C', 'G')])
    num_GU_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('G', 'U'), ('U', 'G')])
    num_AC_pairs = sum(1 for (i, j) in pairs if (S[i], S[j]) in [('A', 'C'), ('C', 'A')])

    num_stacked_quartets = 0
    for (i, j) in pairs:
        if (i + 1, j - 1) in pairs:
            num_stacked_quartets += 1

    num_crossing_pairs = 0
    for (i, j) in pairs:
        for (k, l) in pairs:
            if i < k < j < l:
                num_crossing_pairs += 1

    print("Model 4")
    print("the objective")
    print(m.objVal)
    print("Number of AU pairs")
    print(num_AU_pairs)
    print("Number of GC pairs")
    print(num_GC_pairs)
    print("Number of GU pairs")
    print(num_GU_pairs)
    print("Number of AC pairs")
    print(num_AC_pairs)
    print("Number of Stacked Quartets")
    print(num_stacked_quartets)
    print("Number of Crossing Pairs")
    print(num_crossing_pairs)

    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_plot = 3 * np.cos(theta)
    y_plot = 3 * np.sin(theta)

    x_lab = 4 * np.cos(theta)
    y_lab = 4 * np.sin(theta)

    x_ind = 3.5 * np.cos(theta)
    y_ind = 3.5 * np.sin(theta)

    fig, ax = plt.subplots(figsize=(8, 8))
    circle = plt.Circle((0, 0), 3, fill=False, linewidth=1, linestyle='--')
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-4.5, 4.5)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()
    ax.add_patch(circle)
    ax.scatter(x_plot, y_plot, color='black')

    for i in range(n):
        ax.text(x_lab[i], y_lab[i], S[i], fontsize=12, ha='center', va='center')
        ax.text(x_ind[i], y_ind[i], str(i), ha='center', va='center', fontsize=8, color='red')

    for i, j in pairs:
        ax.plot([x_plot[i], x_plot[j]], [y_plot[i], y_plot[j]], 'g-')

    plt.show()
    return pairs

#Testing!!
S = 'ACGUGCCACGAU'
pairs = solve_nucleic_acid_folding_model_1(S)
pairs = solve_nucleic_acid_folding_model_2(S)
pairs = solve_nucleic_acid_folding_model_3(S)
pairs = solve_nucleic_acid_folding_model_4(S)
print("Pairs formed:", pairs)

S = 'CUUGGCUGGAAACGUAAGUA'
pairs = solve_nucleic_acid_folding_model_1(S)
pairs = solve_nucleic_acid_folding_model_2(S)
pairs = solve_nucleic_acid_folding_model_3(S)
pairs = solve_nucleic_acid_folding_model_4(S)
print("Pairs formed:", pairs)

S = 'GCAUAUGGUCGACGCCUUCAAUAACGAUAC'
pairs = solve_nucleic_acid_folding_model_1(S)
pairs = solve_nucleic_acid_folding_model_2(S)
pairs = solve_nucleic_acid_folding_model_3(S)
pairs = solve_nucleic_acid_folding_model_4(S)
print("Pairs formed:", pairs)

S = 'AGGUACGCCGCUAGAGCGAACCGGGCACCAGCUACGCCGU'
pairs = solve_nucleic_acid_folding_model_1(S)
pairs = solve_nucleic_acid_folding_model_2(S)
pairs = solve_nucleic_acid_folding_model_3(S)
pairs = solve_nucleic_acid_folding_model_4(S)
print("Pairs formed:", pairs)