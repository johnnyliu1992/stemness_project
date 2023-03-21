#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:52:11 2022

@author: Jiannan Liu
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
os.chdir('/Users/johnny/Documents/17p_stemness')

stemness_all=pd.read_csv('cmap_stemness_all.txt', header=0, index_col=0, sep='\t')

inst_info = pd.read_csv("GSE92742_Broad_LINCS_inst_info.txt", sep="\t")
inst_info=inst_info.loc[inst_info['inst_id'].isin(list(stemness_all.index))]
cell_info = pd.read_csv("GSE92742_Broad_LINCS_cell_info.txt", sep="\t")
pert_info=pd.read_csv("GSE92742_Broad_LINCS_pert_info.txt", sep="\t")

#filter out tunor cell lines
cell_info=cell_info.loc[cell_info['sample_type']=='tumor']

#prepare 17p loss/neutral cell lines
cell_17p_loss=pd.read_csv('17p_loss_celllines_all.csv', header=None, index_col=None, sep='\t')
cell_17p_loss=[i.split('_')[0] for i in list(cell_17p_loss.iloc[:,0])]
cell_17p_neutral=pd.read_csv('17p_neutral_celllines_all.csv', header=None, index_col=None, sep='\t')
cell_17p_neutral=[i.split('_')[0] for i in list(cell_17p_neutral.iloc[:,0])]

#find all 17p loss cell lines
celllines_17p_loss=[i for i in list(cell_info['cell_id']) if i in cell_17p_loss]
###### NO 'YAPC'
celllines_17p_loss.remove('YAPC')

#check all instances for cell_17p_loss
inst_info_17p_loss=inst_info.loc[inst_info['cell_id'].isin(celllines_17p_loss)]

#get all compound information
compound_count=inst_info_17p_loss['pert_iname'].value_counts()
compound_list=list(set(inst_info_17p_loss['pert_iname']))

#calculate before treatment stemness for all cell lines
#consider 'DMSO','UnTrt','EMPTY_VECTOR' as un_tre
un_tre_mean={}
for curr_cell in celllines_17p_loss:
    curr_data=inst_info_17p_loss.loc[inst_info_17p_loss['pert_iname'].isin(['DMSO','UnTrt','EMPTY_VECTOR'])]
    curr_data=curr_data.loc[curr_data['cell_id']==curr_cell]
    un_tre_mean[curr_cell]=stemness_all.loc[list(curr_data['inst_id'])].mean().stemness

#init a DataFrame
summary_info=pd.DataFrame(columns=['comp', 'cell', 'inst_num', 'un_tre', 'tre', 'diff', 'fold_change'])

#loop compounds
for curr_comp in compound_list:
    #get all inst recs with curr_comp
    curr_inst=inst_info_17p_loss.loc[inst_info_17p_loss['pert_iname']==curr_comp]
    #loop cell lines
    for curr_cell in list(set(curr_inst['cell_id'])):
        '''
        curr_comp='lacZ'
        curr_cell='HT29'
        '''
        curr_ids=list(curr_inst.loc[curr_inst['cell_id']==curr_cell]['inst_id'])
        curr_mean=stemness_all.loc[curr_ids].mean().stemness
        curr_fold_change=round((curr_mean - un_tre_mean[curr_cell])/un_tre_mean[curr_cell], 3)
        
        summary_info=summary_info.append({'comp':curr_comp, 'cell':curr_cell, 'inst_num':len(curr_ids),
                                          'un_tre':un_tre_mean[curr_cell], 'tre':curr_mean, 
                                          'diff':un_tre_mean[curr_cell]-curr_mean, 
                                          'fold_change':curr_fold_change
                                          }, ignore_index=True)

#for each compound, cell lines >= 5
cummary_comp_count=summary_info['comp'].value_counts()
cummary_comp_count_lg5=cummary_comp_count.loc[cummary_comp_count>=5]
#filtered compound
summary_comp_lg5=summary_info.loc[summary_info['comp'].isin(list(cummary_comp_count_lg5.index))]
#calculate a mean diff for each compound
summary_comp_lg5['comp_mean_diff']=None
comp_mean_diff={i:summary_comp_lg5.loc[summary_comp_lg5['comp']==i]['diff'].mean() for i in list(set(summary_comp_lg5['comp']))}
for i in range(len(summary_comp_lg5)):
    summary_comp_lg5.iloc[i,-1]=comp_mean_diff[summary_comp_lg5.iloc[i,0]]
#sort
summary_comp_lg5=summary_comp_lg5.sort_values(by='comp_mean_diff', ascending=False)

#diff should larger than 0
comp_diff_lg_0=summary_comp_lg5.loc[summary_comp_lg5['diff']>0]
cell_line_counts_comp_final=comp_diff_lg_0['comp'].value_counts()
#find compound insts have more than 5 cell line treated
comp_final=comp_diff_lg_0.loc[comp_diff_lg_0['comp'].isin(list(cell_line_counts_comp_final.loc[cell_line_counts_comp_final>=5].index))]
#calculate min value and sort based on min ==> better plot
comp_final['min']=None
for i in range(len(comp_final)):
    curr_min=comp_final.loc[comp_final['comp']==comp_final.iloc[i,:].comp]['diff'].min()
    comp_final.iloc[i,-1]=curr_min
comp_final=comp_final.sort_values(by='min', ascending=False)
#export data
comp_final.to_csv('conpound_summary_final_loss.csv', index=True, header=True, sep='\t')
