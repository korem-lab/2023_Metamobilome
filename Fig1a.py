from pandas import read_pickle, DataFrame
from mne.stats import fdr_correction
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from utils import read_meta_cnio
from os.path import join

#infiles
meta_file = 'Data/metadata/CNIO_metadata.tsv'
vsgv_file = 'Data/SGVFinder/vsgv_max_spacing_1_del_cooc_0.25msc_75.df'
outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'



def assosciate_variable_region_and_CNIO_data(sample_region, meta ):
    res = {}
    pc_list = meta.loc[meta.subject_disease_status == 'PC'].index.to_list()
    cnt_list = meta.loc[meta.subject_disease_status == 'CTR'].index.to_list()
    for region in sample_region.columns:
        sub = sample_region.dropna(subset=[region])
        temp_pc = sub.loc[sub.index.isin(pc_list), region].to_list()
        temp_ctr = sub.loc[sub.index.isin(cnt_list), region].to_list()
        if len(temp_pc) < 10 or len(temp_ctr) <10:
            continue

        stat, pvalue = mannwhitneyu(temp_pc,temp_ctr, alternative='two-sided')
        effect_size = stat/(len(temp_pc)*len(temp_ctr))
        nums = f'({len(temp_pc)}|{len(temp_ctr)})'
        res[region] = [pvalue, effect_size, nums]

    res_df=DataFrame.from_dict(res, orient='index', columns=['pvalue', 'effect_size', 'numbers'])
    res_df["qvalue"] = fdr_correction(res_df['pvalue'], alpha=0.1, method='indep')[1]
    return res_df

def main():
    #prep metadata
    meta = read_meta_cnio(meta_file)
    meta_2 = meta.loc[meta['subject_disease_status'] != 'Pancreatitis']
    meta_2_stool = meta_2.loc[meta_2['source'] == 'ST']

    #read svfinder output
    vsgv = read_pickle(vsgv_file)

    #analyze disease association
    mw_res_stool = assosciate_variable_region_and_CNIO_data(vsgv, meta_2_stool)
    mw_res_stool[['Bacteria', 'reg_bins']] = mw_res_stool.index.to_series().str.split(':', expand=True)
    mw_res_stool['Bacteria'] = mw_res_stool.Bacteria.replace('445970.PRJNA19655','A. putredinis')
    mw_res_stool['Bacteria'] =mw_res_stool.Bacteria.replace('657318.PRJNA39159','E. rectale')
    mw_res_stool['Bacteria'] = mw_res_stool.Bacteria.replace('1121115.PRJNA195783','B. wexlerae')
    mw_res_stool['Bacteria'] = mw_res_stool.Bacteria.replace('748224.PRJNA46809','F. prausnitzii')
    mw_res_stool['Bacteria'] =  mw_res_stool.Bacteria.replace('411479.PRJNA18195','B. uniformis')
    mw_res_stool['Bacteria'] =mw_res_stool.Bacteria.replace('1235788.PRJNA175977','B. sartorii')
    mw_res_stool['Bacteria'] =mw_res_stool.Bacteria.replace('904306.PRJNA53573','S. vestibularis')
    mw_res_stool['Bacteria'] =mw_res_stool.Bacteria.replace('435591.PRJNA13485','P. distasonis')
    mw_res_stool['Bacteria'] = mw_res_stool.Bacteria.replace('515620.PRJNA29073','E. eligens')

    lines = []
    plt.close('all')
    fig = plt.figure(figsize=(2.5, 2))
    mw_res_stool['trans_pvalue'] = -np.log10(mw_res_stool['pvalue'])
    mw_res_stool.loc[mw_res_stool['qvalue'] >= 0.2, 'Bacteria'] = 'Non Significant'

    col_dict = {}
    val_list = list(sns.color_palette("Paired", 7))
    for key in mw_res_stool.Bacteria.unique():
        for value in val_list:
            col_dict[key] = value
            val_list.remove(value)
            break
    col_dict['Non Significant'] = 'gainsboro'

    g = sns.scatterplot(x='effect_size', y='trans_pvalue', data = mw_res_stool, hue='Bacteria', \
                        palette = col_dict, alpha=0.8, linewidth=0, s=8)

    plt.xlabel('')
    lines.append('xlabel - Effect Size')

    plt.ylabel('')
    lines.append('ylabel - p-value')
    plt.yticks(color='w')
    plt.xlim(0,1)

    plt.hlines(1.553, -12,12, linestyle='dashed', colors='dimgrey', linewidth=0.5) #q 0.2
    plt.hlines(2.38, -12,12, linestyle='dashed', colors='dimgrey', linewidth=0.5) #q 0.1

    g.get_legend().remove()
    g.set(xticklabels=[])
    g.tick_params(bottom=False)

    plt.savefig(join(outdir,'fig1a.pdf'),bbox_inches='tight', dpi=600)



if __name__ == "__main__":
    main()
