import pandas as pd
import sys
import numpy as np
import os
import math
import treeswift as ts

folders = [
    "3435_NT_T1",
    "3513_NT_T1_Fam",
    "3515_Lkb1_T1_Fam",
    "3724_NT_All",
    "3726_NT_T1"
]

data2log = {}
data2branches = {}
data2total_branches = {}
for folder in folders:
    
    target_path = os.path.join('../../result', folder, 'star_cdp_number_of_sol.csv')
    
    tree_path = os.path.join('../../result', folder, 'star_cdp_strict_consensus.tre')
    
    with open(target_path, 'r') as tpr:
        number_of_sols = int(tpr.read())
    
    tre = ts.read_tree_newick(tree_path)
    
    data2branches[folder] = tre.num_nodes(leaves=False, internal=True) - 1
    data2total_branches[folder] = tre.num_nodes(leaves=True, internal=False) - 2

    log_of_sols = math.log(number_of_sols, 2)
    data2log[folder] = log_of_sols


sys.stdout.write("\\begin{table}[!h]\n")
sys.stdout.write("\\caption[")
sys.stdout.write("The logarithm base 2 of the # of optimal tree toplogies]{\\textbf{Tree accuracy and runtime (in seconds) for Startle simulated data sets.} ")
sys.stdout.write("Mean $\\pm$ standard deviations are across replicates for each method. Tree accuracy metrics were computed after SH-contraction of both true and estimated trees.}\n")
sys.stdout.write("\\label{tab:number-of-solutions-shared-branches}\n")
sys.stdout.write("\\centering\n")
sys.stdout.write("\\tiny\n")
sys.stdout.write("\\begin{tabular}{c c c c}\n")
sys.stdout.write("\\toprule \n")
sys.stdout.write("\\DATA of & \\logarithm base 2 of & \\# of branches shared by all trees & \\# of branches in a tree \\\\\n")
sys.stdout.write("\\midrule\n")
for k,v in data2log.items():
    sys.stdout.write(f"\\ {k} & \\ $\\mathbf{{{str(round(v,2))}}}$ & \\ {data2branches[k]} & \\ {data2total_branches[k]}\\\\\n")

sys.stdout.write("\\bottomrule\n")
sys.stdout.write("\\end{tabular}\n")
sys.stdout.write("\\end{table}\n")