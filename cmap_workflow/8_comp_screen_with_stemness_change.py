#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:55:41 2022

@author: Jiannan Liu
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statannot import add_stat_annotation


import os
os.chdir('/Users/johnny/Documents/17p_stemness')

comp_summary_loss=pd.read_csv('conpound_summary_final_loss.csv', index_col=0,header=0,sep='\t')
#change to minus fold change
comp_summary_loss['fold_change']=-comp_summary_loss['fold_change']

comp_summary_loss['std']=None
#calculate std
for i in range(len(comp_summary_loss)):
    comp_summary_loss.iloc[i,-1]=comp_summary_loss.loc[comp_summary_loss['comp']==comp_summary_loss.iloc[i,0]]['diff'].std()
#add a group column
comp_summary_loss['17p Status']='Loss'
    
    
#get all comps
comp_all=list(set(comp_summary_loss['comp']))

#read 17p neutral comp summary
comp_summary_neutral=pd.read_csv('comp_summary_final_neutral.csv', index_col=None,header=0,sep='\t')
#change to minus fold change
comp_summary_neutral['fold_change']=-comp_summary_neutral['fold_change']

comp_summary_neutral=comp_summary_neutral.loc[comp_summary_neutral['comp'].isin(comp_all)]
#add a group column
comp_summary_neutral['17p Status']='Intact'

plot_data_all=pd.concat([comp_summary_loss[['comp','fold_change','17p Status']],comp_summary_neutral[['comp','fold_change','17p Status']]], axis=0)

#calculate t test
ttest_result=pd.DataFrame(columns=['comp','stat','pval'])
for curr_comp in comp_all:
    curr_recs=plot_data_all.loc[plot_data_all['comp']==curr_comp]
    a=list(curr_recs.loc[curr_recs['17p Status']=='Loss']['fold_change'])
    b=list(curr_recs.loc[curr_recs['17p Status']=='Intact']['fold_change'])
    curr_t=stats.ttest_ind(a, b, equal_var=False)
    ttest_result=ttest_result.append({'comp':curr_comp,'stat':curr_t[0],'pval':curr_t[1]}, ignore_index=True)
#filter out by pval
ttest_result=ttest_result.loc[ttest_result['pval']<0.05]
ttest_result=ttest_result.sort_values(by='stat', ascending=False)
#get comp list
comp_list_ttest=list(ttest_result['comp'])[0:35]

#remove some compounds
comp_list_ttest.remove('clonidine')
comp_list_ttest.remove('nimodipine')
comp_list_ttest.remove('RHO-kinase-inhibitor-III[rockout]')



#init filter
plot_data_all=plot_data_all.loc[plot_data_all['comp'].isin(comp_list_ttest)]

#add box pairs
box_pairs=[]   
for curr_comp in comp_list_ttest:
    box_pairs.append(((curr_comp,'Loss'),(curr_comp,'Intact')))

#plot
fig=plt.figure()
#plt.tight_layout()
my_pal={'Loss':'#3953a3', 'Intact':'#EC2025' }
sns.set(rc={'figure.figsize':(20,7), 'axes.facecolor':'black', 'figure.facecolor':'black'})
sns.set_style('white')
#box edge color
flierprops = dict(markeredgecolor='black', markerfacecolor='black')
PROPS = {
    'boxprops':{'edgecolor':'black'},
    'medianprops':{'color':'black'},
    'whiskerprops':{'color':'black'},
    'capprops':{'color':'black'}, 
    'flierprops':flierprops
}


ax = sns.boxplot(x="comp", y="fold_change", hue="17p Status", data=plot_data_all, palette=my_pal, width=0.7, linewidth=1, **PROPS)
ax.set_xlabel('')
ax.set_ylabel('- Fold Change')
ax.set(ylim=(-0.8, 1))
plt.setp(ax.lines, color="0")
plt.xticks(rotation=90)
plt.legend(title='17p Status', labelcolor='black')
#set colors
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['right'].set_color('black')
ax.xaxis.label.set_color('black')
ax.yaxis.label.set_color('black')
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')

#add t test anno
add_stat_annotation(ax, data=plot_data_all, x="comp", y="fold_change", hue="17p Status", box_pairs=box_pairs,
                    test='t-test_welch', loc='inside', verbose=2, comparisons_correction=None, color='black')


#fig.savefig('17p_stemness_35_fc.png', folormat='png', dpi=800, bbox_inches = "tight")
fig.savefig('17p_stemness_32_fc_linewidth_1.pdf', bbox_inches='tight')

#export 35 comps
comps35=pd.DataFrame(list(ttest_result['comp']))
comps35.to_csv('sig_compounds_list.txt', header=False, index=False, sep='\t')

