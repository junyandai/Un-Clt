import pandas as pd
import argparse
import treeswift as ts
import sys
import os
import numpy as np



def set_tree_length(tre):
    for nod in tre.traverse_postorder():
        nod.set_edge_length(1)





def main(args):
    trees1 = ts.read_tree(args.tree1, 'newick')
    trees2 = ts.read_tree(args.tree2, 'newick')
    trees = [trees1] + trees2
    cmat = pd.read_csv(args.cmat, index_col=0, dtype=str)
    taxons = cmat.index.to_list()
    matrix_list = [[] for i in range(len(trees))]
    count = 0
    for tre in trees:
        set_tree_length(tre)
        distance_matrix  = tre.distance_matrix(leaf_labels=True)
        for i in range(len(taxons)):
            for j in range(i + 1, len(taxons)):
                matrix_list[count].append(distance_matrix[taxons[i]][taxons[j]])
        count += 1
    matrix_arr = np.array(matrix_list)
    np.save(args.output, matrix_arr)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t1", "--tree1", type=str,
                        help="File containing reference trees in newick format",
                        required=True)

    parser.add_argument("-t2", "--tree2", type=str,
                        help="File containing other trees in newick format",
                        required=True)

    parser.add_argument('-o', '--output', type=str,
                        help='File containing output result',
                        required=True)

    parser.add_argument('-i', '--cmat', type=str,
                        help='File containing character matrix',
                        required=True)

    main(parser.parse_args())
