import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation
import mygene


cn_data_all=pd.read_csv(r'C:\Users\liuji\Documents\17p_loss\17p_cell_line_ananlysis\CCLE_gene_cn.csv',header=0,sep=',',index_col=0)
cn_data_all_columns=list(cn_data_all.columns)
cn_data_all_columns=[i.split(' ')[0] for i in cn_data_all_columns]
cn_data_all.columns=cn_data_all_columns

#get genes on 17p
genes_17=pd.read_csv(r'C:\Users\liuji\Documents\17p_loss\chromo_17_gene_list.txt',header=0,sep='\t',index_col=0)
genes_17p=[]
for i in range(len(genes_17)):
    if genes_17.iloc[i,1].startswith('17p'):
        genes_17p.append(genes_17.iloc[i,:].name)
genes_17p_final=[i for i in genes_17p if i in cn_data_all_columns]

#get 17p genes' CN data
cn_data_17p=cn_data_all[genes_17p_final]

#read ccle annotation file
ccle_anno=pd.read_csv(r'C:\Users\liuji\Documents\17p_loss\17p_cell_line_ananlysis\ccle_anno.txt',header=0,sep='\t',index_col=1)

#change cell line names of CN data
cn_data_17p=cn_data_17p.loc[list(ccle_anno.index)].dropna()
cn_data_17p_index=list(cn_data_17p.index)
cn_data_17p_index_new=[ccle_anno.loc[i]['CCLE_ID'] for i in cn_data_17p_index]
cn_data_17p.index=cn_data_17p_index_new

all_cells=ccle_anno['CCLE_ID'].to_list()
#ccle_anno.loc[ccle_anno['Site_Primary']=='breast']['CCLE_ID'].to_list()
#get CN data for all on 17p
all_cn_data=cn_data_17p.loc[all_cells].dropna()
all_cn_data['mean']=None
for i in range(len(all_cn_data)):
    all_cn_data.iloc[i,-1]=sum(all_cn_data.iloc[i,0:-1].to_list())/len(all_cn_data.iloc[i,0:-1].to_list())
#get 17p loss cell line list
cell_17p_loss_all=list(all_cn_data.loc[all_cn_data['mean']<0.8].index)
cell_17p_neutral_all=list(all_cn_data.loc[all_cn_data['mean'].between(0.8,1.2)].index)


