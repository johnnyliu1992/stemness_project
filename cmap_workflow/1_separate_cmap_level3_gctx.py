import pandas as pd
from multiprocessing import Pool
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.write_gct as wg

import os
os.chdir('/Users/johnny/Documents/17p_stemness')

def seperate_gct(col_names_list):
    curr_recs=parse('GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx', cid=col_names_list[0:-1])
    wg.write(curr_recs, 'exp_matrix_divided/exp_matrix_'+col_names_list[-1]+'.gct')

'''
This code is used to separate CMAP level3 data into small files to calculate
stemness level

'''

if __name__ == '__main__':
    ''' prepare column ids for gct file  '''
    inst_info = pd.read_csv("GSE92742_Broad_LINCS_inst_info.txt", sep="\t")
    cell_info = pd.read_csv("GSE92742_Broad_LINCS_cell_info.txt", sep="\t")
    
    #find tunro cell lines
    tumor_cells=list(cell_info.loc[cell_info['sample_type']=='tumor']['cell_id'])
    #find column ids based on cell line names
    col_names=list(inst_info.loc[inst_info['cell_id'].isin(tumor_cells)]['inst_id'])
    #divide col_names
    col_names_div=[col_names[i:i+10000] for i in range(0,len(col_names), 10000)]
    #add a index
    for i in range(len(col_names_div)):
        col_names_div[i].append(str(i))
    
    #run multi process job
    with Pool(processes=8) as pool:
        results = pool.map(seperate_gct, col_names_div)
