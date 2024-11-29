import os
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.text import Text
import seaborn as sns

proj_dir = '../../'
datas = ['3435_NT_T1','3513_NT_T1_Fam', '3515_Lkb1_T1_Fam', '3724_NT_All', '3726_NT_T1']

dfs = []
for data in datas:
    cur_mig_df = pd.read_csv(os.path.join(proj_dir, 'result', data, 'RF_distances.csv'), sep=',')
    cur_mig_df['data'] = data
    dfs.append(cur_mig_df)

all_df = pd.concat(dfs, ignore_index=True)

print(all_df)

fig = plt.figure(figsize=(14.4, 8))
sns.boxplot(x='data', \
                y='fprate', \
                data=all_df)
plt.savefig('fprate.pdf', dpi=300, bbox_inches='tight', format='pdf')
