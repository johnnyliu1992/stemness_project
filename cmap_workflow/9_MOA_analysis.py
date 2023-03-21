#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 20:59:27 2022

@author: Jiannan Liu
"""
import pandas as pd
import os
os.chdir('/Users/johnny/Documents/17p_stemness')

'''  MoA  '''
moa_data=pd.read_csv('Repurposing_Hub_export.txt', index_col=0,header=0,sep='\t')
#comp_summary=pd.read_csv('sig_compounds_list.txt', index_col=None,header=None,sep='\t')
comp_summary=pd.read_csv('35_compounds_list.txt', index_col=None,header=None,sep='\t')
comp_list=list(comp_summary.iloc[:,0])
#filter out not in moa_data
comp_list_filtered=[i for i in comp_list if i in list(moa_data.index)]
#init df
moa_sum=pd.DataFrame(columns=['comp','moa'])
for curr_comp in comp_list_filtered:
    curr_rec=moa_data.loc[curr_comp]
    if ', ' in curr_rec['MOA']:
        curr_moa_list=curr_rec['MOA'].split(', ')
        for i in curr_moa_list:
            moa_sum=moa_sum.append({'comp':curr_comp,'moa':i}, ignore_index=True)
    else:
        moa_sum=moa_sum.append({'comp':curr_comp,'moa':curr_rec['MOA']}, ignore_index=True)

moa_count=pd.DataFrame(moa_sum['moa'].value_counts())
moa_count['comp']=None
for i in range(len(moa_count)):
    moa_count.iloc[i,-1]=', '.join(moa_sum.loc[moa_sum['moa']==moa_count.iloc[i,:].name]['comp'].to_list())
moa_count=moa_count.reset_index()
moa_count.columns=['moa','comp_count','comp']

#moa_count.to_csv('112_compound_moa_summary.txt', index=False, header=True, sep='\t')
moa_sum.to_csv('25_compound_moa_summary.tsv', index=False, header=True, sep='\t')
