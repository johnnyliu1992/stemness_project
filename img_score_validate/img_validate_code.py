import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import seaborn as sns
from scipy import stats
from statannot import add_stat_annotation
import matplotlib.pyplot as plt



import os
os.chdir('/Users/johnny/Documents/17p_stemness/img_score_validate')


img_score=pd.read_csv('Full_Rev_TNBC_is_norm_False_epoch_100.csv', index_col=None, header=0, sep=',')
img_score=img_score.drop_duplicates(subset=['pid'], keep='first')
img_score=img_score.sort_values(by='pred_score')


brca_exp=pd.read_csv('brca_mrna', index_col=0, header=0, sep='\t')
brca_exp=brca_exp[[i for i in list(brca_exp.columns) if i.endswith('01')]]
brca_exp.columns=[i[8:12] for i in list(brca_exp.columns)]

all_pids=list(img_score['pid'])

final_exp_quart=brca_exp[all_pids[0:20] + all_pids[-20::]]
final_exp_quart.to_csv('img_score_vali_quart_exp_gsea.txt',index=True, header=True, sep='\t')


final_exp_med=brca_exp[all_pids[0:40] + all_pids[-40::]]
final_exp_med.to_csv('img_score_vali_med_exp_gsea.txt',index=True, header=True, sep='\t')





'''

stemness maker correlation

'''
img_score=pd.read_csv('Full_Rev_TNBC_is_norm_False_epoch_100.csv', index_col=None, header=0, sep=',')
img_score=img_score.drop_duplicates(subset=['pid'], keep='first')

markers=['WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','PLCB1','ITGA6','GPR161','MYC','PRKACA','MAPK8','HIF1A','CD44','EZH2','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','DHH','PTCH1','PRKCA']
#init a df
stemness_corr=pd.DataFrame(columns=['gene','corr','pval'])
all_pids=list(img_score['pid'])
all_pre_scores=list(img_score['pred_score'])

for curr_marker in markers:
    curr_exp=list(brca_exp[all_pids].loc[curr_marker])
    curr_corr=spearmanr(all_pre_scores,curr_exp)
    stemness_corr=stemness_corr.append({'gene':curr_marker,'corr':curr_corr[0],'pval':curr_corr[1]}, ignore_index=True)


'''

stemness maker boxplot
'''
img_score=pd.read_csv('Full_Rev_TNBC_is_norm_False_epoch_100.csv', index_col=None, header=0, sep=',')
img_score=img_score.drop_duplicates(subset=['pid'], keep='first')
img_score=img_score.sort_values(by='pred_score')
all_pids=list(img_score['pid'])

pids_l=all_pids[0:20]
pids_h=all_pids[-20::]

markers=['WNT4','MDM2','PROM1','HHAT','CAMK2B','KLF4','NOTCH1','EPAS1','GJA1','PLCB1','ITGA6','ABCG2','GPR161','MYC','PRKACA','MAPK8','BMI1','DVL1','CTNNB1','SOX9','HIF1A','CD24','CD44','EZH2','RUNX2','TWIST1','PPP3CA','ROCK2','HHIP','LRP6','BCL2','EPCAM','KDM5B','FSTL1','NFATC1','NES','CDON','CD34','DHH','PTCH1','PRKCA']
for curr_marker in markers:
    curr_exp_l=list(brca_exp[pids_l].loc[curr_marker])
    curr_exp_h=list(brca_exp[pids_h].loc[curr_marker])
    curr_plot_data=pd.DataFrame()
    curr_plot_data['low']=curr_exp_l
    curr_plot_data['high']=curr_exp_h
    
    #plot violin plot
    fig=plt.figure()
    sns.set(rc={'figure.figsize':(3,5)})
    sns.set_style('white')
    ax = sns.boxplot(data=curr_plot_data, width=0.7, linewidth=2)
    ax.set_xlabel('')
    ax.set_title(curr_marker)
    
    #add t test anno
    add_stat_annotation(ax, data=curr_plot_data,box_pairs=[('low','high')],
                        test='t-test_welch', loc='inside', verbose=2, comparisons_correction=None)
    
    #fig.savefig('cmap_17p_stemness_before_treat.png', format='png', dpi=800, bbox_inches = "tight")
    fig.savefig('plots/'+curr_marker+'.pdf')    
    


'''
prepare gene rank list for fgsea (R package)

'''
rank_data=pd.read_csv('ranked_gene_list_stemness_markers.tsv', index_col=0, header=0, sep='\t')
rank_data=rank_data[['SCORE']]
#rank_data.to_csv('ranked_gene_list_stemness_markers.rnk', header=False, index=False, sep='\t')


id_map=pd.read_csv('gene_id_map.txt', index_col=0, header=0, sep='\t')
id_map=id_map[~id_map.index.duplicated(keep='first')]
rank_data=pd.concat([rank_data,id_map], axis=1)
rank_data=rank_data[['ENTREZID','SCORE']]
rank_data=rank_data.dropna()
rank_data.to_csv('ranked_gene_list_stemness_markers_entrez.rnk', header=True, index=False, sep='\t')



'''
investigate CN profile using PMID: 35705804 data

'''
cn_profile=pd.read_csv('cnv_profile_tcga.txt', index_col=1, header=0, sep='\t')
cn_profile_brca=cn_profile.loc[cn_profile['cancer_type']=='BRCA']
cn_profile_brca.index=[i[-4::] for i in list(cn_profile_brca.index)]


img_score=pd.read_csv('Full_Rev_TNBC_is_norm_False_epoch_100.csv', index_col=None, header=0, sep=',')
img_score=img_score.drop_duplicates(subset=['pid'], keep='first')
img_score=img_score.sort_values(by='pred_score')


all_pids=list(img_score['pid'])
all_pids.remove('A1EW')
all_pids.remove('AA10')

cn_profile_hpci=cn_profile_brca.loc[all_pids[0:40]+all_pids[-40::]]

data_l=list(cn_profile_hpci['ploidy'].loc[all_pids[0:40]])
data_h=list(cn_profile_hpci['ploidy'].loc[all_pids[-40::]])
curr_plot_data=pd.DataFrame()
curr_plot_data['HPCI low']=data_l
curr_plot_data['HPCI high']=data_h

#plot violin plot

fig=plt.figure()
#plt.tight_layout()
my_pal={'HPCI low':'#3953a3', 'HPCI high':'#EC2025' }
sns.set(rc={'figure.figsize':(3,4)})
sns.set_style('white')
ax = sns.boxplot(data=curr_plot_data, palette=my_pal, width=0.7, linewidth=2)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('ploidy')

#add t test anno
add_stat_annotation(ax, data=curr_plot_data,box_pairs=[('HPCI low','HPCI high')],
                    test='t-test_welch', loc='inside', verbose=2, comparisons_correction=None)
fig.savefig('LOH.pdf')


'name', 'cancer_type', 'sex', 'barcodeTumour', 'barcodeNormal',
   'tumour_mapd', 'normal_mapd', 'GC_correction_before',
   'GC_correction_after', 'RT_correction_before', 'RT_correction_after',
   'n_het_SNP', 'n_segs_logR', 'n_segs_BAF', 'n_segs_logRBAF_diff',
   'frac_homo', 'purity', 'ploidy', 'goodness_of_fit', 'n_segs',
   'segs_size', 'n_segs_1kSNP', 'homdel_segs', 'homdel_largest',
   'homdel_size', 'homdel_fraction', 'LOH', 'mode_minA', 'mode_majA',
   'WGD', 'GI', 'QC'



'''

compare updated HPCI score

'''

img_score=pd.read_csv('Full_Rev_TNBC_is_norm_False_epoch_100.csv', index_col=None, header=0, sep=',')
img_score=img_score.drop_duplicates(subset=['pid'], keep='first')
img_score=img_score.sort_values(by='pred_score')

img_score_l=list(img_score['pid'])[0:40]
img_score_h=list(img_score['pid'])[-40::]


img_score_new=pd.read_csv('OCT4_20_0.45.csv', index_col=None, header=0, sep=',')
img_score_new.columns=['pid', 'TNBC_label', 'TNBC_score', 'TNBC_predicted_score']
img_score_new=img_score_new.drop_duplicates(subset=['pid'], keep='first')
img_score_new=img_score_new.sort_values(by='TNBC_predicted_score')

img_score_new_l=list(img_score_new['pid'])[0:40]
img_score_new_h=list(img_score_new['pid'])[-40::]


ol_l=[i for i in img_score_l if i in img_score_new_l]
ol_h=[i for i in img_score_h if i in img_score_new_h]

















