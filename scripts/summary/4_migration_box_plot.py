import os
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.text import Text
import seaborn as sns



proj_dir = '../../'
datas = ['3513_NT_T1_Fam', '3515_Lkb1_T1_Fam', '3724_NT_All']

dfs = []
for data in datas:
    cur_mig_df = pd.read_csv(os.path.join(proj_dir, 'result', data, 'migration_table.csv'), sep=',')
    cur_mig_df['data'] = data
    dfs.append(cur_mig_df)

all_df = pd.concat(dfs, ignore_index=True)

print(all_df)

fig = plt.figure(figsize=(8, 6))
sns.boxplot(x='data', \
                y='migrations', \
                data=all_df)
plt.savefig('migrations.pdf', dpi=300, bbox_inches='tight', format='pdf')

fig = plt.figure(figsize=(8, 6))
sns.boxplot(x='data', \
                y='reseedings', \
                data=all_df)
plt.savefig('reseedings.pdf', dpi=300, bbox_inches='tight', format='pdf')