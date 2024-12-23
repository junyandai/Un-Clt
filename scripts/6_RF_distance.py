import argparse
import pandas
import treeswift
import os
import sys



def infer_states_under_star_at_internal_helper(states_below):
    # Assume ambiguous state = -1
    #        unedited state = 0
    #        edited states = 1, 2, ...

    states_below_list = sorted(list(states_below))

    nbelow = len(states_below_list)
    if nbelow == 1:
        state = states_below_list[0]
    elif nbelow == 2:
        if -1 in states_below:
            state = states_below_list[1]
        else:
            state = 0
    else:
        state = 0

    return state


def infer_states_under_star_at_leaf(node, cmat):
    
    states_at_leaf = cmat.loc[node.label].values.tolist()
    # states_at_leaf = cmat.loc[int(node.label)].values.tolist()
    
    node.states = []
    node.states_below = []
    for state in states_at_leaf:
        state_set = set([state])
        node.states_below.append(state_set)
        node.states += list(state_set)


def infer_states_under_star_at_internal(node, children):
    nchar = len(children[0].states_below)
    node.states = []
    node.states_below = []

    for i in range(nchar):
        # Get states below
        states_below = children[0].states_below[i]
        for child in children[1:]:
            states_below = states_below.union(child.states_below[i])
        node.states_below.append(states_below)

        # Infer ancestral states
        state = infer_states_under_star_at_internal_helper(states_below)
        node.states += [state]

    for child in children:
        del child.states_below


def get_mutations_on_edge(head, tail):
    nchar = len(head.states)
    muts = []
    for i in range(nchar):
        hs = head.states[i]
        ts = tail.states[i]
        if ts == -1:
            pass
        elif hs == ts:
            pass
        else:
            if hs != 0:
                sys.exit("Error in star labeling!")
                muts.append((i, ts))
    return muts


def get_mutations_above_root(root):
    nchar = len(root.states)
    muts = [] 
    for i in range(nchar):
        rs = root.states[i]
        if rs == -1:
            pass
        elif rs == 0:
            pass
        else:
            muts.append((i, rs))
    return muts


def infer_muts_under_star_model(tree, cmat):
    for node in tree.traverse_postorder():
        if node.is_leaf():
            infer_states_under_star_at_leaf(node, cmat)
        else:
            # Check if node is binary
            children = node.child_nodes()

            infer_states_under_star_at_internal(node, children)

            # Annotate child edges with mutations
            for child in children:
                child.muts = get_mutations_on_edge(node, child)
                child.set_edge_length(len(child.muts))

    # Lastly process the root!
    root = tree.root
    root.muts = get_mutations_above_root(root)
    root.set_edge_length(len(root.muts))

    return tree


def contract_zero_length_branches(tree):
    # Set incoming edges incident to root and leaves
    # to be greater than 0 so they are not contracted
    rlen = tree.root.get_edge_length()
    if rlen == 0:
        tree.root.set_edge_length(0.5)

    for node in tree.traverse_leaves():
        if node.get_edge_length() == 0:
            node.set_edge_length(0.5)

    # Contract zero length edges
    for node in tree.traverse_postorder():
        if node.get_edge_length() == 0:
            node.contract()

    # Re-set edge lengths associated with root and leaves
    for node in tree.traverse_leaves():
        if node.get_edge_length() == 0.5:
            node.set_edge_length(0.0)

    if rlen == 0:
        tree.root.set_edge_length(0.0)



# def check_rooted_at_fake(tree, fake):
#     root = tree.root
#     children = root.child_nodes()
#     for child in children:
#         if child.is_leaf():
#             if child.label == fake:
#                 return True
#     return False


def read_character_data(inputf):
    chars = pandas.read_csv(inputf)

    colmap = {}
    for col in chars.columns.values.tolist():
        if col == "Unnamed: 0":
            colmap[col] = "cell"
        else:
            colmap[col] = col
    chars.rename(columns=colmap, inplace=True)

    return chars


def get_shared_leaves(tree1, tree2, fake=None):
    lset1 = set([x for x in tree1.labels()])
    lset2 = set([x for x in tree2.labels()])
    keep = list(lset1.intersection(lset2))
    return keep


def get_clade_strings_at_child_nodes(node, root=None):
    clades = []

    if node.is_leaf():
        # No clades below leaf so return empty list
        return clades

    children = node.child_nodes()
    if len(children) == 1:
        sys.exit("Failed to suppress unifurcations!\n")

    for child in children:
        leaves = [leaf.label for leaf in child.traverse_leaves()]

        if len(leaves) == 1:
            # Found trivial clade so do nothing
            pass
        else:
            leaf_set = set(leaves)
            if root is not None:
                if root in leaf_set:
                    # Found clade that to exclude with re-rooting so do nothing
                    pass
                else:
                    clade = ",".join(sorted(list(leaf_set)))
                    clades.append(clade)
            else:
                clade = ",".join(sorted(list(leaf_set)))
                clades.append(clade)

    return clades


def get_nontrivial_clade_strings(tree, root=None):
    all_clades = []
    for node in tree.traverse_postorder():
        if node.is_leaf():
            # No children below leaf
            pass
        else:
            clades = get_clade_strings_at_child_nodes(node, root=root)
            if len(clades) > 0:
                all_clades += clades
    return all_clades


def compare_trees(tree1, tree2, root=None):
    """
    Compares non-trivial clade set for trees

    If 'root' taxon is specified, then trees are re-rooted at leaf;
    otherwise trees are treated as being rooted correctly.

    Nothing about this function or the helpers is smart but whatever

    Parameters
    ----------
    tree1 : treeswift tree object
            First tree (typically the model tree)
    tree2 : treeswift tree object
            Second tree (typically the estimated tree)
    root : label of leaf on which to root tree

    Returns
    -------
    nl : int
         Size of the shared leaf set, i.e., the number of leaves in both trees
    i1 : int
          Number of non-trivial clades in tree 1 after restriction to shared leaves
    i2 : int
          Number of non-trivial clades in tree 2 after restriction to shared leaves
    fn : int
         Number of non-trivial clades in tree 1 that are not in tree 2
    fp : int
         Number of non-trivial clades in tree 2 that are not in tree 1

    Example
    -------
    If tree 1 corresponds to "(((A,B,C),D),E);"
    and tree 2 corresponds to "((((A,B),C),D),E);",
    In this example, if the root is None or E, then
      + tree 1 and tree 2 share five leaves (A, B, C, D, E).
      + tree 1 has two non-trivial clades "A,B,C" and "A,B,C,D"
      + tree 2 has three non-trivial clades "A,B", "A,B,C", "A,B,C,D"
      + no non-trivial clades in tree 1 are missing from tree 2
      + one non-trivial clade in tree 2 is missing from the tree 1
    so the output is "5 2 3 0 1".

    If the root is A, then 
    
    """
    nl = len([leaf for leaf in tree1.traverse_leaves()])

    clades1 = get_nontrivial_clade_strings(tree1, root=root)
    clades2 = get_nontrivial_clade_strings(tree2, root=root)
    clades2 = set(clades2)

    i1 = len(clades1)
    i2 = len(clades2)

    clade_set2 = set(clades2)

    tp = 0
    for c1 in clades1:
        if c1 in clade_set2:
            tp += 1

    fn = i1 - tp
    fp = i2 - tp

    return (nl, i1, i2, fn, fp, tp)


def main(args):


    table = {'pair':[], 'fn':[], 'fp':[], 'tp':[], 'fnrate':[], 'fprate':[],'tprate':[]}
    trees1 = treeswift.read_tree(args.tree1, 'newick')

    # Read tree 2 and contract edges as specified
    trees2 = treeswift.read_tree(args.tree2, 'newick')
    #print(tree2)

    trees = [trees1] + trees2
    count = 0
    for i in range(len(trees)):
        for j in range(i+1, len(trees)):
            tree1 = trees[i]
            tree2 = trees[j]

            # Restrict trees to their shared leaf set
            keep = get_shared_leaves(tree1, tree2)
            tree1 = tree1.extract_tree_with(keep, suppress_unifurcations=True)
            tree2 = tree2.extract_tree_with(keep, suppress_unifurcations=True)

            #Check trees have root
            if args.root is not None:
                if args.root not in keep:
                    sys.exit("At least one tree is not rooted at %s!" % args.root)

            # Compare trees
            [nl, i1, i2, fn, fp, tp] = compare_trees(tree1, tree2, root=args.root) 
            fnrate = float(fn) / i1 if i1 != 0 else float('inf')
            fprate = float(fp) / i2 if i2 != 0 else float('inf')
            tprate = float(tp) / i2 if i2 != 0 else float('inf')
            table['pair'].append(f'({i},{j})')
            table['fn'].append(fn)
            table['fp'].append(fp)
            table['tp'].append(tp)
            table['fnrate'].append(fnrate)
            table['fprate'].append(fprate)
            table['tprate'].append(tprate)
            print(f'round {count} : {fn},{fp},{tp},')
            count += 1

    print(table)
    table_df = pandas.DataFrame(table)
    table_df.to_csv(args.output,index=False)

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

    # parser.add_argument("-c1", "--contract1", type=int, default=0,
    #                     help="Method for contracting branches in tree 1:\n\n" 
    #                             + "  0 = do not contract (default)\n"
    #                             + "  1 = contract branches with no mutations, using mutations inferred under star homoplasy model (requires -m option)\n"
    #                             + "  2 = contract branches with no mutations, using mutations from input newick string")
    # parser.add_argument("-ex1", "--extract_tree1", type=int, default=0, help="0 = do nothing\n" + "1 = extract newick string for tree 1 from nwk to tre\n");
    
    # parser.add_argument("-ex2", "--extract_tree2", type=int, default=0, help="0 = do nothing\n" + "1 = extract newick string for tree 2 from nwk to tre\n");
    
    # parser.add_argument("-c2", "--contract2", type=int, default=0,
    #                     help="Method for contracting edges in tree 2")

    # parser.add_argument("-m", "--chars", type=str, default=None,
    #                     help="Character data in CSV format")

    parser.add_argument("-r", "--root", type=str, default=None,
                        help="Name of *leaf* node representing unedited root cell")

    # parser.add_argument('-t2p', '--tree2_path', type=str, default=None, help='File to write the extract tree 2')

    main(parser.parse_args())