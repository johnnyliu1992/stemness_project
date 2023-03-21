import pandas as pd
from gtfparse import read_gtf

import os
os.chdir('/Users/johnny/Documents/17p_stemness/rnaseq_analysis')

'''
combine result files
'''
file_list=os.listdir('data')
data_combined=pd.DataFrame(index=list(pd.read_csv('data/'+file_list[0], index_col=0, header=0, sep='\t', skiprows=1).index))

for curr_file in file_list:
    curr_data=pd.read_csv('data/'+curr_file, index_col=0, header=0, sep='\t', skiprows=1)
    curr_sample=curr_data.columns[-1].split('/')[-2]
    curr_data=pd.DataFrame(curr_data.iloc[:,-1])
    curr_data.columns=[curr_sample]
    data_combined=pd.concat([data_combined,curr_data], axis=1)

data_combined.to_csv('rnaseq_result_combined_probe.csv', index=True, header=True, sep='\t')

'''
transform ids to gene symbol
'''
ref=read_gtf("gencode.v36.annotation.gtf")
ref=ref.set_index('gene_id')
#part of ref
ref=ref['gene_name']
ref=ref.loc[~ref.index.duplicated(keep='first')]
ref=pd.DataFrame(ref)
data_combined=pd.concat([data_combined,ref],axis=1)
data_combined=data_combined[list(data_combined.columns)[::-1]]

'''
process duplicated probes
'''
#get all duplicated rows
dup_genes_count=data_combined['gene_name'].value_counts()
dup_genes=list(dup_genes_count.loc[dup_genes_count>=2].index)
data_dup=data_combined.loc[data_combined['gene_name'].isin(dup_genes)]
data_none_dup=data_combined.loc[~data_combined['gene_name'].isin(dup_genes)]

data_dup_new=pd.DataFrame(columns=list(data_dup.columns))
for curr_gene in dup_genes:
    curr_rec=data_dup.loc[data_dup['gene_name']==curr_gene]
    curr_rec=curr_rec.sum()
    curr_rec['gene_name']=curr_gene
    curr_rec=pd.DataFrame(curr_rec).T
    data_dup_new=pd.concat([data_dup_new,curr_rec],axis=0)

#combine with data_none_dup
data_final=pd.concat([data_none_dup, data_dup_new], axis=0)
data_final=data_final[['gene_name','TH_1_S92_L002', 'TH_2_S93_L002', 'TH_3_S94_L002', 'TH_4_S95_L002', 'TH_5_S96_L002', 'TH_6_S97_L002', 'TH_7_S98_L002', 'TH_8_S99_L002', 'TH_9_S100_L002', 'TH_10_S101_L002', 'TH_11_S102_L002', 'TH_12_S103_L002']]
data_final.to_csv('rnaseq_result_combined_symbol.csv', index=False, header=True, sep='\t')


'''
process data for GSEA

'''
#=====>>>> before treat
before_treat_result=pd.read_csv('DE_results_BeforeTreat.txt', index_col=0, header=0, sep='\t')
before_treat_result_sig=before_treat_result.loc[before_treat_result['adj.P.Val']<0.01]
before_treat_result_sig=pd.DataFrame(before_treat_result_sig['t'])
before_treat_result_sig.to_csv('GSEA_data_BeforeTreat.rnk', index=True, header=False, sep='\t')

#for DAVID
before_treat_result=pd.read_csv('DE_results_BeforeTreat.txt', index_col=0, header=0, sep='\t')
before_treat_result_sig=before_treat_result.loc[before_treat_result['adj.P.Val']<0.01]
#before_treat_result_sig=before_treat_result_sig.loc[~before_treat_result_sig['logFC'].between(-1.2,1.2)]
before_treat_result_sig_list_up=pd.DataFrame(list(before_treat_result_sig.loc[before_treat_result_sig['t']>0].index))
before_treat_result_sig_list_up.to_csv('david_results/data/de_gene_list_BeforeTreat_up.txt', index=False, header=False, sep='\t')
before_treat_result_sig_list_down=pd.DataFrame(list(before_treat_result_sig.loc[before_treat_result_sig['t']<0].index))
before_treat_result_sig_list_down.to_csv('david_results/data/de_gene_list_BeforeTreat_down.txt', index=False, header=False, sep='\t')


#=====>>>> 17p treat
loss17_treat_result=pd.read_csv('DE_results_Loss17pTreat.txt', index_col=0, header=0, sep='\t')
loss17_treat_result_sig=loss17_treat_result.loc[loss17_treat_result['adj.P.Val']<0.05]
loss17_treat_result_sig=pd.DataFrame(loss17_treat_result_sig['t'])
loss17_treat_result_sig.to_csv('GSEA_data_Loss17pTreat.rnk', index=True, header=False, sep='\t')

#for DAVID
loss17_treat_result=pd.read_csv('DE_results_Loss17pTreat.txt', index_col=0, header=0, sep='\t')
loss17_treat_result_sig=loss17_treat_result.loc[loss17_treat_result['adj.P.Val']<0.05]
loss17_treat_result_sig=loss17_treat_result_sig.loc[~loss17_treat_result_sig['logFC'].between(-1.2,1.2)]
loss17_treat_result_sig_list_up=pd.DataFrame(list(loss17_treat_result_sig.loc[loss17_treat_result_sig['t']>0].index))
loss17_treat_result_sig_list_up.to_csv('david_results/data/de_gene_list_17ploss_treat_up.txt', index=False, header=False, sep='\t')
loss17_treat_result_sig_list_down=pd.DataFrame(list(loss17_treat_result_sig.loc[loss17_treat_result_sig['t']<0].index))
loss17_treat_result_sig_list_down.to_csv('david_results/data/de_gene_list_17ploss_treat_down.txt', index=False, header=False, sep='\t')

#=====>>>> diff
diff_result=pd.read_csv('DE_results_Diff.txt', index_col=0, header=0, sep='\t')
diff_result_sig=diff_result.loc[diff_result['adj.P.Val']<0.05]
diff_result_sig=pd.DataFrame(diff_result_sig['t'])
diff_result_sig.to_csv('GSEA_data_Diff.rnk', index=True, header=False, sep='\t')

#for DAVID
diff_result=pd.read_csv('DE_results_Diff.txt', index_col=0, header=0, sep='\t')
diff_result_sig=diff_result.loc[diff_result['adj.P.Val']<0.05]
diff_result_sig=diff_result_sig.loc[~diff_result_sig['logFC'].between(-1.2,1.2)]
diff_result_sig_list_up=pd.DataFrame(list(diff_result_sig.loc[diff_result_sig['t']>0].index))
diff_result_sig_list_up.to_csv('david_results/data/de_gene_list_diff_up.txt', index=False, header=False, sep='\t')
diff_result_sig_list_down=pd.DataFrame(list(diff_result_sig.loc[diff_result_sig['t']<0].index))
diff_result_sig_list_down.to_csv('david_results/data/de_gene_list_diff_down.txt', index=False, header=False, sep='\t')



###=====>>>> top 1500 DE up and down  ====>> before treat
before_treat_result=pd.read_csv('DE_results_BeforeTreat.txt', index_col=0, header=0, sep='\t')
before_treat_result_sig=before_treat_result.loc[before_treat_result['adj.P.Val']<0.05]
before_treat_result_sig=before_treat_result_sig.sort_values(by='logFC', ascending=False)

before_treat_result_sig_list_up=pd.DataFrame(list(before_treat_result_sig.index)[0:2000])
before_treat_result_sig_list_up.to_csv('david_results/data/de_gene_list_BeforeTreat_up.txt', index=False, header=False, sep='\t')
before_treat_result_sig_list_down=pd.DataFrame(list(before_treat_result_sig.index)[-2000::])
before_treat_result_sig_list_down.to_csv('david_results/data/de_gene_list_BeforeTreat_down.txt', index=False, header=False, sep='\t')


###=====>>>> top 1500 DE up and down  ====>> 17p loss
loss17_treat_result=pd.read_csv('DE_results_Loss17pTreat.txt', index_col=0, header=0, sep='\t')
loss17_treat_result_sig=loss17_treat_result.loc[loss17_treat_result['adj.P.Val']<0.05]

loss17_treat_result_sig=loss17_treat_result_sig.sort_values(by='logFC', ascending=False)
loss17_treat_result_sig_list_up=pd.DataFrame(list(loss17_treat_result_sig.index)[0:2000])
loss17_treat_result_sig_list_up.to_csv('david_results/data/de_gene_list_17ploss_treat_up.txt', index=False, header=False, sep='\t')
loss17_treat_result_sig_list_down=pd.DataFrame(list(loss17_treat_result_sig.index)[-2000::])
loss17_treat_result_sig_list_down.to_csv('david_results/data/de_gene_list_17ploss_treat_down.txt', index=False, header=False, sep='\t')


###=====>>>> top 1500 DE up and down  ====>> diff
diff_result=pd.read_csv('DE_results_Diff.txt', index_col=0, header=0, sep='\t')
diff_result_sig=diff_result.loc[diff_result['adj.P.Val']<0.05]

diff_result_sig=diff_result_sig.sort_values(by='logFC', ascending=False)
diff_result_sig_list_up=pd.DataFrame(list(diff_result_sig.index)[0:2000])
diff_result_sig_list_up.to_csv('david_results/data/de_gene_list_diff_up.txt', index=False, header=False, sep='\t')
diff_result_sig_list_down_recs=diff_result_sig.loc[diff_result_sig['logFC']<0]
diff_result_sig_list_down=pd.DataFrame(list(diff_result_sig_list_down_recs.index))
diff_result_sig_list_down.to_csv('david_results/data/de_gene_list_diff_down.txt', index=False, header=False, sep='\t')


stemness_markers=['CD44','CD24','EPCAM','RUNX2','SOX2','SOX9','ABCG2','EZH2','ITGA6','MYC','GJA1','ALDH1A1','FSTL1','LRP6','ROCK2','MAPK8','FZD10','PLCB1','CAMK2A','PPP3CA','PRKCA','NFATC1','DVL1','WNT4','MDM2','CAMK2B','HHAT','DHH','GPR161','PTCH1','HHIP','BCL2','PRKACA','CDON','TWIST1','NOTCH1','KLF4','CD44','ABCG2','CD34','EZH2','HIF1A','PROM1','MYC','EPAS1','BMI1','KDM5B','NES','CTNNB1']
overlapped=[i for i in stemness_markers if i in list(diff_result.index)]
#heatmap_rec=diff_result.loc[overlapped].sort_values(by='adj.P.Val', ascending=True)
heatmap_rec=diff_result.loc[overlapped].sort_values(by='logFC', ascending=False)
heatmap_rec=heatmap_rec.drop_duplicates(keep='first')
genes_remove=['GJA1','ABCG2','BMI1','DVL1','CTNNB1','SOX9','CD24','RUNX2','TWIST1','PPP3CA','CD34']
heatmap_rec=heatmap_rec.loc[~heatmap_rec.index.isin(genes_remove)]
overlapped_sorted=list(heatmap_rec.index)




#volvano plot prepare
before_treat_result=pd.read_csv('DE_results_BeforeTreat.txt', index_col=0, header=0, sep='\t')
before_treat_result_sig=before_treat_result.loc[before_treat_result['adj.P.Val']<0.05]
stemness_markers=['WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','PLCB1','ITGA6','GPR161','MYC','PRKACA','MAPK8','HIF1A','CD44','EZH2','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','DHH','PTCH1','PRKCA']
sig_genes=[i for i in stemness_markers if i in list(before_treat_result_sig.index)]
sig_rec=before_treat_result_sig.loc[sig_genes].sort_values(by='logFC', ascending=False)
sig_rec_list=list(sig_rec.index)

#volvano plot data mod
before_treat_result=pd.read_csv('DE_results_BeforeTreat.txt', index_col=0, header=0, sep='\t')
stemness_markers=['WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','PLCB1','ITGA6','GPR161','MYC','PRKACA','MAPK8','HIF1A','CD44','EZH2','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','DHH','PTCH1','PRKCA']
nan_markers=['BST1','NT5E','NMRK1','QPRT','ENPP3','SIRT1']
mod_index=list(before_treat_result.index)+stemness_markers+nan_markers
before_treat_result_mod=before_treat_result.loc[mod_index]
before_treat_result_mod=before_treat_result_mod.drop_duplicates(keep='last')
before_treat_result_mod.to_csv('DE_results_BeforeTreat_volcano.txt', header=True, index=True, sep='\t')


'''
David result plot
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os
os.chdir('/Users/johnny/Documents/17p_stemness/rnaseq_analysis')

plot_data=pd.read_csv('david_plot_data_diff_down.csv', index_col=None, header=0, sep=',')
plot_data['Term']=[i.split(':')[1] for i in list(plot_data['Term'])]
plot_data['Fold Enrichment']=[round(i,2) for i in list(plot_data['Fold Enrichment'])]
plot_data=plot_data.sort_values(by='Fold Enrichment', ascending=False)

fig=plt.figure()

#matplotlib.rc_file_defaults()
sns.set_style('white')
fig, ax1 = plt.subplots(figsize=(7,5))
sns.barplot(data = plot_data, 
            x='Term', 
            y='Fold Enrichment', 
            color="#c81c26",
            alpha=0.9,
            ax=ax1)
ax1.set_ylim([0, 5])
#set colors
ax1.spines['bottom'].set_color('black')
ax1.spines['top'].set_color('black')
ax1.spines['left'].set_color('black')
ax1.spines['right'].set_color('black')
ax1.xaxis.label.set_color('black')
ax1.yaxis.label.set_color('black')
ax1.tick_params(axis='x', colors='black')
ax1.tick_params(axis='y', colors='black')


ax2 = ax1.twinx()
ax1.tick_params(axis='x', rotation=90)
ax1.tick_params(axis='y', rotation=90)
ax1.bar_label(ax1.containers[0], color='black')

sns.lineplot(data=plot_data, 
             x='Term', 
             y='PValue', 
             marker='o', 
             sort=False, 
             ax=ax2)
ax2.set_ylim([0, 0.1])
ax2.tick_params(axis='y', rotation=90)
#set colors
ax2.spines['bottom'].set_color('black')
ax2.spines['top'].set_color('black')
ax2.spines['left'].set_color('black')
ax2.spines['right'].set_color('black')
ax2.xaxis.label.set_color('black')
ax2.yaxis.label.set_color('black')
ax2.tick_params(axis='x', colors='black')
ax2.tick_params(axis='y', colors='black')

#fig.savefig('david_plot_diff_down.png', format='png', dpi=800, bbox_inches = "tight")
fig.savefig('david_plot_diff_down.pdf', bbox_inches='tight')





#plot
fig=plt.figure()
sns.set(rc={'figure.figsize':(3,5)})
sns.set_style('white')
ax=sns.barplot(data=plot_data, x='Count', y='Pathways', color="#c81c26", orient = 'h')
ax.bar_label(ax.containers[0])
fig.savefig('MOA_summary_plot.png', format='png', dpi=800, bbox_inches = "tight")






'''
DE data process, chaneg order of genes

'''
soure_data=pd.read_csv('DE_results_BeforeTreat_volcano.txt', index_col=0, header=0, sep='\t')

l_3marker=['DHH','CAMK2B','WNT4','BCL2','CDON','HHAT','PTCH1','PRKACA','HHIP','MDM2','EPAS1','GPR161','MAPK8','FSTL1','KDM5B','NOTCH1','LRP6','PROM1', 'OCT4','SOX2','NANOG']















