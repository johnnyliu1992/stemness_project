''' id and format transfer ''' 

import pandas as pd
from multiprocessing import Pool
import os
os.chdir('/Users/johnny/Documents/17p_stemness')

def filetrans(params):#params[file_list,symbol_list]
    for curr_file in params[0]:
        data=pd.read_csv('exp_matrix_divided/'+curr_file, header=0, index_col=0, sep='\t', skiprows=2)
        data.index=params[1]
        data.to_csv('exp_matrix_divided_symbol/'+curr_file, index=True, header=True, sep='\t')



if __name__ == '__main__':    
    #read one dataset
    data = pd.read_csv('exp_matrix_divided/exp_matrix_0.gct', header=0, index_col=0, sep='\t', skiprows=2)
    #read two ref files
    gene_ref=pd.read_csv('gene_id_ref.txt', header=0, index_col=0, sep='\t', converters={'entrezgene_id':str})
    gene_ref_1=pd.read_csv('id_transfer_ref.csv', header=0, index_col=None, sep='\t', converters={'entrez':str})
    
    gene_list=list(data.index)
    #transfer entrez to gene symbol
    gene_list_symbol=[]
    for curr_gene in gene_list:
        if len(gene_ref.loc[gene_ref['entrezgene_id']==str(curr_gene)]['hgnc_symbol'])!=0:
            gene_list_symbol.append(list(gene_ref.loc[gene_ref['entrezgene_id']==str(curr_gene)]['hgnc_symbol'])[0])
        else:
            try:
                gene_list_symbol.append(list(gene_ref_1.loc[gene_ref_1['entrez']==str(curr_gene)]['symbol'])[0])
            except:
                gene_list_symbol.append(str(curr_gene))
                #print(curr_gene)
    #assemble a params
    all_files=os.listdir('/Users/johnny/Documents/17p_stemness/exp_matrix_divided')
    exist_files=os.listdir('/Users/johnny/Documents/17p_stemness/exp_matrix_divided_symbol')
    all_files=[i for i in all_files if i not in exist_files]
    all_files.remove('.DS_Store')
    
    all_files_div=[all_files[i:i+1] for i in range(0,len(all_files), 1)]
    params=[[i,gene_list_symbol] for i in all_files_div]
    #run multi process job
    with Pool(processes=7) as pool:
        results = pool.map(filetrans, params)