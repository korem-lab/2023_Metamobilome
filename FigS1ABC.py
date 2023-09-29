from pandas import read_pickle, DataFrame
from mne.stats import fdr_correction
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from utils import add_p_val, read_meta_cnio
from skbio.stats.ordination import pcoa
from skbio.stats.distance import anosim, permanova, DistanceMatrix
from scipy.spatial.distance import squareform, pdist
from os.path import join

#infiles
meta_file = 'Data/metadata/CNIO_metadata.tsv'
vsgv_file = 'Data/SGVFinder/vsgv_max_spacing_1_del_cooc_0.25msc_75.df'
outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def main():

#prep metadata
    meta = read_meta_cnio(meta_file)
    meta_2 = meta.loc[meta['subject_disease_status'] != 'Pancreatitis']
    meta_2_stool = meta_2.loc[meta_2['source'] == 'ST']

    #read svfinder output
    vsgv = read_pickle(vsgv_file)
    vsgv2= vsgv.loc[vsgv.index.isin(meta_2_stool.index),].T
    vsgv2['genome'] = vsgv2.index.to_series().str.split(':').str[0]
    res =[]
    variable='subject_disease_status'
    beta_div_metric = 'canberra'
    #make pcoa figures for the two example genomes
    for genome in vsgv2.genome.unique():
        if len(vsgv2.loc[vsgv2['genome'] == genome].dropna(axis=1, how='all').columns) - 1 > 9:
            ordf = vsgv2.loc[vsgv2['genome'] == genome].dropna(axis=1, how='all').drop(columns=['genome']).T
            temp_meta = ordf.merge(meta_2_stool, left_index=True, right_index=True, how='left')
            mat = squareform(pdist(ordf.to_numpy(), metric = beta_div_metric))

            permanova_res = permanova(DistanceMatrix(mat), temp_meta[variable], permutations=999)
            text_str = f'for {genome} permanova pval = {str(permanova_res[5])}'
            res.append([genome, ordf.shape[0], permanova_res[5]])
            if genome in ['1121115.PRJNA195783', '445970.PRJNA19655']:
                xa, ya = 'PC1', 'PC2'
                plt.close('all')
                pcoa_mat = pcoa(mat)
                df = pcoa_mat.samples.copy()
                df.index = ordf.index
                x2 = df.merge(temp_meta, left_index=True, right_index=True, how = 'left')

                plt.figure(figsize=(2,2))
                g = sns.scatterplot(x=xa,y=ya, data = x2, hue = variable, s=15, palette=dict(PC="salmon", CTR="darkcyan"))
                plt.xlabel("")
                plt.ylabel("")
                g.set(xticklabels=[])
                g.set(yticklabels=[])
                g.get_legend().remove()
                print(f'{genome} perm pval {str(permanova_res[5])}')
                print('%s (%.2f%%)' % (xa, pcoa_mat.proportion_explained.loc[xa]*100))
                print('%s (%.2f%%)' % (ya, pcoa_mat.proportion_explained.loc[ya]*100))
                plt.savefig(f'{outdir}{genome}_vsv_pcoa.pdf', dpi=300)

    res_df = DataFrame(res, columns=['Bacteria', 'num_samples', 'permanova p-val'])
    res_df['permanova q-val'] = fdr_correction(res_df['permanova p-val'], alpha=0.1, method='indep')[1]
    res_df['Bacteria'] = res_df.Bacteria.replace('445970.PRJNA19655','A. putredinis')
    res_df['Bacteria'] =res_df.Bacteria.replace('657318.PRJNA39159','E. rectale')
    res_df['Bacteria'] = res_df.Bacteria.replace('1121115.PRJNA195783','B. wexlerae')
    res_df['Bacteria'] = res_df.Bacteria.replace('748224.PRJNA46809','F. prausnitzii')
    res_df['Bacteria'] =  res_df.Bacteria.replace('411479.PRJNA18195','B. uniformis')
    res_df['Bacteria'] =res_df.Bacteria.replace('1235788.PRJNA175977','B. sartorii')
    res_df['Bacteria'] =res_df.Bacteria.replace('904306.PRJNA53573','S. vestibularis')
    res_df['Bacteria'] =res_df.Bacteria.replace('435591.PRJNA13485','P. distasonis')
    res_df['Bacteria'] = res_df.Bacteria.replace('515620.PRJNA29073','E. eligens')
    res_df['Bacteria'] = res_df.Bacteria.replace('1128111.PRJNA80237','V. atypica')
    res_df['Bacteria'] = res_df.Bacteria.replace('1316408.PRJNA196120','Streptococcus sp. HSISM1')
    res_df['Bacteria'] = res_df.Bacteria.replace('879309.PRJNA52053','Veillonella sp. oral taxon 158 str. F0412')
    res_df.drop(columns=['num_samples']).to_csv(join(outdir,'FigS1A_svs_indiv_bacteria_permanova_res.csv'))





if __name__ == "__main__":
    main()
