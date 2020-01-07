""" trim down the epi cells in wrkdir/scVCF_filtered_all/ to just the 
	tumor cells, ie. the only ones it actually needs to be run on """
import os

import pandas as pd

vcf_dir = '/home/ubuntu/cerebra/cerebra/wrkdir/scVCF_filtered_all/'
tumor_cells_df = pd.read_csv('tumor_cells_list.csv', header=None,
                             names=['name'])
tumor_cells_list = list(tumor_cells_df['name'])

for f in os.listdir(vcf_dir):
    cell = f.strip('.vcf')
    if cell not in tumor_cells_list:
        os.remove(vcf_dir + f)
# print(cell)
