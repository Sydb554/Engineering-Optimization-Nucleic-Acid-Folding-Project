# ISyE 3133 Project: The Nucleic Acid Folding Problem
## Overview
The Nucleic Acid Folding Problem aims to predict the secondary structure of a nucleic acid molecule given only its nucleotide sequence. For this project, we represent the nucleotide sequence as a circular string. The secondary structure is determined by predicting which nucleotides will pair together within the circle, where each nucleotide can only form one pair. An example of a nucleotide pairing is illustrated below:

Figure 1: Pairing {(0, 11),(1, 4),(3, 10),(2, 6)}


## Project Description
1. The First Crude Model
The first model assumes:

Only nested pairings (non-crossing pairs) are allowed
Only complementary pairs (A, U) and (G, C) are considered
The most stable pairing is the one with the most matched pairs
Figure 2: A nested pairing with complementary nucleotides


Objective: Formulate an integer linear program to find the most likely pairing under these assumptions.

2. Simple Biological Enhancements
In the second model, we enhance the biological assumptions:

Only nested pairings are allowed
Any pair of nucleotides must be at least three positions away from each other
Complementary pairs (A, U) and (G, C) are considered
Some non-complementary pairs (G, U) and (A, C) are allowed
The most stable pairing is the most likely
The stability of a pairing is a weighted function of each pair type
Table 1: Weights of different pairs in Model 2

Pair	Weight
(G, C)	3
(A, U)	2
(G, U)	0.1
(A, C)	0.05
Objective: Formulate an integer linear program to find the most likely pairing under these assumptions.

3. More Complex Biological Enhancements
In the third model, we include stacked quartets:

Any pair of nucleotides must be at least three positions away from each other
Complementary pairs (A, U) and (G, C) are considered
Some non-complementary pairs (G, U) and (A, C) are allowed
The stability of a pairing is a weighted function of each pair type and the number of stacked quartets
Table 2: Weights of pairs and stacked quartets in Model 3

Item	Weight
(G, C) pair	3
(A, U) pair	2
(G, U) pair	0.1
(A, C) pair	0.05
Stacked quartet	1
Objective: Formulate an integer linear program to find the most likely pairing under these assumptions.

4. Model with Crossing Pairs
The most complex model allows for crossing pairs. Here, we adjust Model 3 to allow up to 10 crossing pairs.

Figure 3: Example of crossing pairs


## Report
The final report includes:

Introduction: An overview of the problem and background information.
Model: Explanation of the modeling approach, including variables, objectives, and constraints.
Solution and Analysis: Solutions for each model using the provided dataset, analyzing running time and model complexity.
Final Project Details
Note: This project assumes no prior knowledge of biology, only understanding of integer linear programming formulations and modular/clock arithmetic.

Final project for Engineering Optimization (ISyE 3133), completed Spring 2024.

## Code Structure
models/: Contains the linear programming models for each problem variation.
data/: Contains the input nucleotide strings.
scripts/: Python scripts for data processing and model execution.
images/: Visual representations of nucleotide pairing and pairing examples.

## Requirements
Python 3.8+
Gurobi (or any other linear programming solver)
Numpy
Matplotlib
