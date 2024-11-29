import pandas as pd 
import numpy as np
import os
import sys


def binarize(df, priors=None):
    columns = list(df.columns)
    cell_list = df[columns[0]].values
    new_rows = {}
    new_cols = [columns[0]]
    for cell in cell_list:
        new_rows[cell] = [cell]
    weights = []

    for col in columns[1:]:
        character = df[col].values

        state_set = set(character)
        try:
            state_set.remove(0)
        except KeyError:
            pass
        try:
            state_set.remove(-1)
        except KeyError:
            pass

        states = sorted(list(state_set))
        states = np.array(states)
        nstates = states.size

        for state in states:
            xdf = priors[(priors["character"] == col) & \
                         (priors["state"] == state)]

            if xdf.shape[0] != 1:
                sys.exit("Bad mutation priors!\n")

            prob = xdf.probability.values[0]
            wght = -1.0 * np.log(prob)
            wght = int(np.round(np.round(wght, 3), 2) * 100)
            weights.append(str(wght))

        new_cols += [col + '_' + str(state) for state in states]

        for cell, x in zip(cell_list, character):
            if x == 0:
                new_rows[cell] += [0] * nstates
            elif x == -1:
                new_rows[cell] += [-1] * nstates
            else:
                data = [0] * nstates
                indx = np.where(states == x)[0][0]
                data[indx] = 1
                new_rows[cell] += data

    cols = new_cols
    rows = [row for row in new_rows.values()]
    return pd.DataFrame(rows, columns=cols), weights


def write_to_nexus(df:pd.DataFrame, paup_file:str, folder: str, binary_file:str, taxon_map_file:str, weights:list):
    columns = list(df.columns)
    nchar = len(columns) - 1
    ntax = df.shape[0]

    with open(binary_file, 'w') as fp1, \
         open(taxon_map_file, 'w') as fp2:
        fp1.write("#NEXUS\n\n")
        fp1.write("Begin data;\n")
        fp1.write("\tDimensions ntax=%d nchar=%d;\n" % (ntax + 2, nchar))
        fp1.write("\tFormat datatype=standard gap=-;\n")
        fp1.write("\tMatrix\n")

        fake = ''.join(['0'] * nchar)
        fp1.write("ROOT0\n")
        fp1.write(fake + '\n')
        fp1.write("ROOT1\n")
        fp1.write(fake + '\n')

        for index, row in df.iterrows():
            cell = row.values[0]
            new_label = str("LEAF%d" % index)
            fp1.write(new_label + '\n')
            fp2.write("%s,%s\n" % (cell, new_label))

            data = row.values[1:]
            data = [str(x) for x in data]
            data = ''.join(data)
            data = data.replace('-1', '-')
            fp1.write(data + '\n')

        fp1.write('\t;\n')
        fp1.write('End;\n')

    # Note - 
    # irreversible (Camin-Sokal); up means higher numbers are derived
    # you can use weight set of integers...
    # I can probably do round to second decimal and take base 100 weights
    all_trees_file = os.path.join(folder, 'all_trees.trees')
    with open(paup_file, 'w') as fp:
        fp.write("#NEXUS\n")
        fp.write("BEGIN PAUP;\n")
        fp.write("set autoclose=yes warntree=no warnreset=no;\n")
        fp.write(f"execute {binary_file};\n")
        fp.write("outgroup ROOT0;\n")
        fp.write("typeset myctype = irrev.up:1-%d;\n" % nchar)
        fp.write("wtset mywtset vector = " + ' '.join(weights) + ";\n")
        fp.write("assume typeset=myctype wtset=mywtset;\n")
        fp.write("hsearch start=stepwise addSeq=random swap=None nreps=10 rseed=55555;")
        fp.write("hsearch start=1 swap=TBR nbest=500 rseed=12345;\n")
        fp.write("rootTrees;\n")
        fp.write(f"savetrees File={all_trees_file} root=yes trees=all format=newick;\n")
        fp.write("END;\n")


def main():
    
    proj_dir = "/fs/cbcb-lab/ekmolloy/jdai123/un-clt-study/"

    data_dir = os.path.join(proj_dir, 'data')

    result_dir = os.path.join(proj_dir, 'result')
    
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
        
        priors = pd.read_csv(priors_path, sep=',')
        cmat_df = pd.read_csv(cmat_path, sep = ',')
        binary_file = os.path.join(cur_res_path, 'binary.nex')

        taxon_map_file = os.path.join(cur_res_path, 'taxon_map.csv')
        paup_file = os.path.join(cur_res_path, 'paup_camsok_hsearch_fast.nex')
        # print(cmat_df)
        cmat_df.astype({col: int for col in cmat_df.columns[1:]})
        binary_df, weights = binarize(cmat_df, priors)
        write_to_nexus(binary_df, paup_file, cur_res_path, binary_file, taxon_map_file, weights)
        print("Finished writing paup nex file for " + cur_res_path)
        



if __name__ == '__main__':
    main()