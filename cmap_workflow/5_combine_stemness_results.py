'''  stemness ressult combine   '''  

import pandas as pd
import os
os.chdir('/Users/johnny/Documents/17p_stemness')


all_files=os.listdir('/Users/johnny/Documents/17p_stemness/final_stemness')
all_files.remove('.DS_Store')

all_stemness=pd.DataFrame(columns=['sample', 'stemness'])
for curr_file in all_files:
    curr_data=pd.read_csv('final_stemness/'+curr_file, header=None, index_col=None, sep='\t')
    curr_data.columns=['sample', 'stemness']
    all_stemness=pd.concat([all_stemness, curr_data], ignore_index=True)

#export combined final data
all_stemness.to_csv('cmap_stemness_all.txt', index=False, header=True, sep='\t')
