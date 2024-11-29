import os
import subprocess as sp


def write_star_cdp_sbatch(exe: str, cmat: str, priors: str, sbatch_file: str, search_space: str):
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
                    ' '.join([time,'-v','-o','star_cdp_usage.log',exe,'-i',cmat, '-m', priors, '-x 0','-t',search_space, '-consensus','-rand 1000','-XOUTG','--output', 'star_cdp',\
                    '&>','star_cdp.log','2>&1'])
                    ]

    with open(sbatch_file, 'w') as sf:
        sf.write("\n".join(rows_to_wrtie))


def submit_star_cdp_sbacth(sbatch: str, data_prefix: str, res_dir:str):
    job_name = "--job-name=" + "c.star_cdp." + data_prefix
    output =  "--output="+ "c.star_cdp." + data_prefix + ".%j.out"
    err = "--error=" + "c.star_cdp." + data_prefix + ".%j.err"
    os.chdir(res_dir)
    if not os.path.exists(os.path.join(res_dir,'star_cdp_one_sol.tre')):
        sp.run(['sbatch', job_name, output, err, sbatch])

def main():
    
    proj_dir = '/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study/'
    result_dir = os.path.join(proj_dir, 'result')

    data_dir = os.path.join(proj_dir, 'data')
    
    software_dir = os.path.join(proj_dir, 'software')

    star_cdp_exe = os.path.join(software_dir, 'Star-CDP','src','star-cdp')

    folders = [f for f in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, f))]

    

    for folder in folders:
        cur_data_path = os.path.join(data_dir, folder)
        cur_res_path = os.path.join(result_dir, folder)
        
        if not os.path.exists(cur_res_path):
            os.mkdir(cur_res_path)
        
            
        cmat_path = os.path.join(cur_res_path, folder + '_pruned_character_matrix_with_outg.csv')
        priors_path = os.path.join(cur_data_path, folder + '_priors.csv')
        print(priors_path)
            
        search_space = os.path.join(cur_res_path, 'search_space.trees')
        star_cdp_sbatch = os.path.join(cur_res_path, 'run_star_cdp.sbatch')
        data_prefix = folder + "/"
        print(data_prefix)

            #need to change remind! star_cdp_one_sol.tre
        if not os.path.exists(os.path.join(cur_res_path, "star_cdp_one_sol.tre")):
            write_star_cdp_sbatch(star_cdp_exe, cmat_path, priors_path, star_cdp_sbatch, search_space)
            print('Writing sbatch file for ' + cur_res_path)
            submit_star_cdp_sbacth(star_cdp_sbatch, data_prefix, cur_res_path)
            print('Submiting  job for ' + cur_res_path)
        else:
            print(f"already exists for {folder}, so gonna skip it")
                
            

if __name__ == '__main__':

    main()