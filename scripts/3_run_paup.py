import os
import subprocess as sp




def write_paup_sbatch(curr_dir : str, paup_exe: str, paup_nex_file:str, sbatch_file:str):
    # sbatch_file_name = 'run_nni.sbatch'
    # sbatch_file = os.path.join(curr_dir, sbatch_file_name)
    time="/usr/bin/time"
    rows_to_wrtie = ['#!/bin/bash', 
                     '#SBATCH --time=24:00:00', 
                     '#SBATCH --cpus-per-task=1', 
                     '#SBATCH --ntasks=1',
                    '#SBATCH --mem=48G',
                    '#SBATCH --qos=highmem',
                    '#SBATCH --partition=cbcb',
                    '#SBATCH --account=cbcb',
                    '#SBATCH --constraint=EPYC-7313',
                    '#SBATCH --exclusive',
                    'module load Python3/3.8.15',
                    ' '.join([time, '-v', '-o', 'paup_usage.log', paup_exe, paup_nex_file, '&>', 'paup.log', '2>&1'])
                    ]

    with open(sbatch_file, 'w') as sf:
        sf.write("\n".join(rows_to_wrtie))



def submit_paup_sbacth(sbatch: str, data_prefix: str, res_dir:str):
    job_name = "--job-name=" + "b.paup." + data_prefix
    output =  "--output="+ "b.paup." + data_prefix + ".%j.out"
    err = "--error=" + "b.paup." + data_prefix + ".%j.err"
    os.chdir(res_dir)
    if not os.path.exists(os.path.join(res_dir, 'all_trees.trees')):
        sp.run(['sbatch', job_name, output, err, sbatch])

def main():
    
    proj_dir = '/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study/'
    result_dir = os.path.join(proj_dir, 'result')

    data_dir = os.path.join(proj_dir, 'data')
    
    software_dir = os.path.join(proj_dir, 'software')

    paup_exe = os.path.join(software_dir, 'paup4a168_centos64')

    folders = [f for f in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, f))]

    for folder in folders:
        cur_data_path = os.path.join(data_dir, folder)
        cur_res_path = os.path.join(result_dir, folder)
        
        if not os.path.exists(cur_res_path):
            os.mkdir(cur_res_path)
        

            
        cmat_path = os.path.join(cur_data_path, folder + '_pruned_character_matrix.csv')
        priors_path = os.path.join(cur_data_path, folder + '_priors.csv')   
        print(priors_path)
            
        paup_nex_file = os.path.join(cur_res_path, 'paup_camsok_hsearch_fast.nex')
        paup_sbatch = os.path.join(cur_res_path, 'run_paup.sbatch')
        data_prefix = folder + "/"
        print(data_prefix)
            

        if not os.path.exists(os.path.join(cur_res_path, "all_trees.trees")):
            write_paup_sbatch(cur_res_path, paup_exe, paup_nex_file, paup_sbatch)
            print('Writing sbatch file for ' + cur_res_path)
            submit_paup_sbacth(paup_sbatch, data_prefix, cur_res_path)
            print('Submiting job for ' + cur_res_path)
                
            

if __name__ == '__main__':

    main()