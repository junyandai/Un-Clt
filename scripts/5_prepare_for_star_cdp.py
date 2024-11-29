import os
import subprocess as sp
import shutil
import pandas as pd


def add_outgroup(newick, outgroup='0'):
    """
    Adds an outgroup to a Newick string.
    The Newick tree is assumed to end with a semicolon.
    The outgroup is added by modifying the base call with a new root that includes the outgroup.
    """
    # Remove the last semicolon and then add the outgroup
    newick_trimmed = newick.strip()[:-1]  # Remove semicolon and any extra whitespace
    return f"({outgroup},{newick_trimmed});"

def process_newick_file(input_file, output_file, outgroup):
    """
    Reads a file with Newick strings, adds an outgroup to each, and writes them to a new file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            updated_newick = add_outgroup(line, outgroup)
            outfile.write(updated_newick + '\n')



def main():
    proj_dir = '/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study/' 
    data_path = os.path.join(proj_dir, 'data')
    
    result_dir = os.path.join(proj_dir, 'result')
    
    paup_dir = os.path.join(proj_dir, 'result')

    folders = [f for f in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, f))]
    
   
    for folder in folders:
        cur_paup_path = os.path.join(paup_dir, folder)
        cur_res_path = os.path.join(result_dir, folder)
        cur_data_path = os.path.join(data_path, folder)

        if not os.path.exists(cur_res_path):
            os.mkdir(cur_res_path)

            
        paup_usage_path = os.path.join(cur_paup_path, 'paup_usage.log')

        paup_search_space_path = os.path.join(cur_paup_path, 'paup_trees.trees')

        outgroup_camt_path = os.path.join(cur_res_path, folder + '_pruned_character_matrix_with_outg.csv')

        if not os.path.exists(outgroup_camt_path) or True:
            cmat_path = os.path.join(cur_data_path, folder + '_pruned_character_matrix.csv')
            cmat_df = pd.read_csv(cmat_path, index_col=[0])
            out_cmat_df = pd.DataFrame([{col: 0 for col in cmat_df.columns}], columns=cmat_df.columns[0:], index=[0])
            out_cmat_df = pd.concat([out_cmat_df, cmat_df], ignore_index=False)
            out_cmat_df.to_csv(outgroup_camt_path)


        data_prefix = folder + "/"
        print(data_prefix)

        if not os.path.exists(os.path.join(cur_res_path, 'search_space.trees')):
            process_newick_file(paup_search_space_path, os.path.join(cur_res_path, 'search_space.trees'), '0')    
        
        if not os.path.exists(os.path.join(cur_res_path, 'paup_usage.log')):
            shutil.copy(paup_usage_path, os.path.join(cur_res_path, 'paup_usage.log'))

        


        
                
            

if __name__ == '__main__':

    main()