import argparse
import pandas as pd
import treeswift as ts
import os
import sys




def count_migration(tre, labels_list, taxon2order):
    optimal_labeling = {}
    score = {}
    trace_back = {}
    count = 0
    internalnode2label = {}
    first_primary_tumor_id = min(labels_list)
    for node in tre.traverse_postorder():
        score[node] = {}
        trace_back[node] = {}
        if not node.is_leaf():
            if not node.is_root():
                internalnode2label[node] = f'I{count}'
                count += 1
                node.label = internalnode2label[node]
            else:
                node.label = 'root'
                internalnode2label[node] = 'root'

        for lab in labels_list:
            if node.is_leaf():
                if str(node.get_label()) not in taxon2order:
                    print("Warning Missing some taxon in dict!!!")
                if lab == taxon2order[str(node.get_label())]:
                    score[node][lab] = 0
                else:
                    score[node][lab] = float('inf')
            else:
                tot_cost = 0
                trace_back[node][lab] = {}

                for child in node.child_nodes():
                    #if len(node.child_nodes()) != 2:
                        #print(f"Warring non-binary tree!! len: {len(node.child_nodes())}")


                    best_child_lab = []
                    best_child_contr = float('inf')
                    
                    for child_lab in labels_list:

                        cost = 1 if child_lab != lab else 0
                        
                        if score[child][child_lab] == float('inf'):
                            continue

                        if best_child_contr > cost + score[child][child_lab] and score[child][child_lab] != float('inf'):
                            best_child_contr = cost + score[child][child_lab]
                            best_child_lab = [child_lab]
                        elif best_child_contr == cost + score[child][child_lab] and best_child_contr != float('inf'):
                            best_child_lab.append(child_lab)

                    tot_cost += best_child_contr
                    trace_back[node][lab][child] = best_child_lab
                
                score[node][lab] = tot_cost
        #print(f'{score[node]}||{node.label}')
    # print(f'The number of internal nodes(include root): {count-1}')
    # print(f'The number of total nodes: {len(list(tre.traverse_postorder()))}')
    return score[tre.root][first_primary_tumor_id],trace_back


def top_down(tre,trace_back,labels_list,order2site, site2order):
    # print(len(list(tre.traverse_postorder())))
    # print(len(trace_back))
    
    for node in tre.traverse_postorder():
        if node not in trace_back:
            print(f"Node label: {node.label} missing in trace_back")
    que = []
    
    first_primary_tumor_id = min(labels_list)
    
    que.append((tre.root,order2site[first_primary_tumor_id]))
    
    optimal_label = {tre.root:order2site[first_primary_tumor_id]}
    
    migration_count = {'p2n_transition':0, 'n2n_transition' : 0, 'n2p_transition':0, 'not_transition':0, 'p2p_transition':0, 'migration':0, 'reseeding':0}
    
    migration_count_all = {}
    primary = order2site[first_primary_tumor_id]
    while que:
        node,lab = que.pop(0)
        for child in node.child_nodes():
            # print(trace_back[node][child])
            
            optimal_label[child] = order2site[min(trace_back[node][site2order[lab]][child])]
            # print(optimal_label[child])
            if lab != optimal_label[child]:
                if f'{lab}->{optimal_label[child]}' not in migration_count_all:
                    migration_count_all[f'{lab}->{optimal_label[child]}'] = 1
                    
                else:
                    migration_count_all[f'{lab}->{optimal_label[child]}'] += 1
                    
                
                # if lab[:-2] != optimal_label[child][:-2] and f'{lab[:-2]}->{optimal_label[child][:-2]}' not in migration_count_all:
                #     migration_count_all[f'{lab[:-2]}->{optimal_label[child][:-2]}'] = 1
                    
                # elif lab[:-2] != optimal_label[child][:-2] and f'{lab[:-2]}->{optimal_label[child][:-2]}' in migration_count_all:
                #     migration_count_all[f'{lab[:-2]}->{(optimal_label[child])[:-2]}'] += 1


                if lab == primary and optimal_label[child] != primary:
                    migration_count['p2n_transition'] += 1
                elif lab != primary and optimal_label[child] != primary:
                    migration_count['n2n_transition'] += 1
                elif lab != primary and optimal_label[child] == primary:
                    migration_count['n2p_transition'] += 1
                    migration_count['reseeding'] += 1
                else:
                    print(f"warnning: parent lab: {lab}; child lab: {optimal_label[child]}")
                migration_count['migration'] += 1
                
            else:
                migration_count['not_transition'] += 1
                # print(f'parent lab: {lab} == child lab {optimal_label[child]}')
            
            que.append((child, optimal_label[child]))
    node_label2site_label = {}
    for node,lab in optimal_label.items():
        node_label2site_label[node.label] = lab
    return node_label2site_label,migration_count,migration_count_all




def main(args):
    table = {'id':[], 'migrations':[],'reseedings':[]}
    trees1 = ts.read_tree(args.tree1, 'newick')

    # Read tree 2 and contract edges as specified
    trees2 = ts.read_tree(args.tree2, 'newick')
    #print(tree2)

    df_kp_meta = pd.read_csv('KPTracer_meta.csv', index_col = 0)

    tumor_table_for_dp = pd.read_csv('tumor_info_table_for_dp.csv', sep=',')

    cmat = pd.read_csv(args.cmat, index_col=0, dtype=str)
    taxon2site = df_kp_meta.loc[cmat.index]['SubTumor']

    for cell, site in taxon2site.items():
        if len(site.split('_')) > 3:
            taxon2site[cell] = '_'.join(site.split('_')[:-1])
    
    cur_info_table_for_dp = tumor_table_for_dp[tumor_table_for_dp['family'] == args.data]

    index_for_dp_list = cur_info_table_for_dp['index_for_dp'].tolist()


    # KPtracer paper define site label -> human readable site label
    name_to_full_site_name_map = dict(zip(cur_info_table_for_dp['name'], cur_info_table_for_dp['full_site_name']))

    # human readable site label -> index 
    name_to_index_map = dict(zip(cur_info_table_for_dp['full_site_name'], cur_info_table_for_dp['index_for_dp']))
    #
    index_to_site_map = dict(zip(cur_info_table_for_dp['index_for_dp'], cur_info_table_for_dp['full_site_name']))
    
    # leaf cell name -> human readable site label 
    for cell,site in taxon2site.items():
        taxon2site[cell] = name_to_full_site_name_map[site]

    taxon2order = {}
    for cell,site in taxon2site.items():
        taxon2order[cell] = name_to_index_map[site]


    all_sites_list = cur_info_table_for_dp['full_site_name'].to_list()


    trees = [trees1] + trees2
    count = 0
    for tre in trees:
        table['id'].append(count)
        count += 1
        score,trace_back = count_migration(tre, index_for_dp_list, taxon2order)
        node2site_label, migration_count, migration_count_all = top_down(tre,trace_back,index_for_dp_list,index_to_site_map, name_to_index_map)
        table['migrations'].append(migration_count['migration'])
        table['reseedings'].append(migration_count['reseeding'])

    result_df = pd.DataFrame(table)
    result_df.to_csv(args.output, index=False)

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

    parser.add_argument('-d', '--data', type=str,
                        help='File containing character matrix',
                        required=True)

    main(parser.parse_args())