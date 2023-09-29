import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from mne.stats import fdr_correction
from pandas import read_csv, read_pickle, DataFrame, concat
from os.path import join

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def get_correl_with_host(signif_df, mge_cnt, mts_log):
    all_res = DataFrame()
    mge_cnt = mge_cnt.set_index('sample_id')
    for cur_mge in signif_df.OTU:
        temp =mts_log.T.merge(mge_cnt[cur_mge].T, left_index=True, right_index=True, how='inner')
        res = []
        for col in temp.columns:
            r, p = spearmanr(temp[cur_mge], temp[col])
            res.append([cur_mge, col, r, p])
        res_df = DataFrame(res[0:-1], columns = ['mge', 'bac', 'spear_r', 'p'])
        res_df['q_val'] = fdr_correction(res_df['p'], alpha=0.1, method='indep')[1]
        all_res = concat([all_res, res_df])
    return all_res


def read_meta_cnio(meta_file):
    meta = read_csv(meta_file, sep='\t')
    meta['source'] = meta.sample_alias.str.split('-').str[2]
    meta = meta.set_index('sample_alias')
    meta = meta[~meta.index.duplicated(keep='first')]
    meta = meta.loc[meta.subject_disease_status != 'NonElegible']
    meta = meta.dropna(subset=['diabin'])
    meta.diabin = meta.diabin.astype('int')
    return meta

def main():

    mts = read_csv('Data/mOTUs/CNIO_motus2.tsv', sep='\t')
    mts = mts.loc[:, mts.columns != 'MMPC-3109233-ST-0']
    s_mts = mts.set_index('#consensus_taxonomy')
    s_mts = s_mts.fillna(0)
    s_mts.drop('-1', inplace=True)
    s_mts = s_mts.truediv(other = s_mts.sum(0), axis = 1)
    s_mts = s_mts.loc[(s_mts != 0).sum(1) > 5]
    mts_log = np.log10(s_mts)
    mts_log = mts_log.replace(-np.Inf, -6)

    plas_cnts = read_csv('Data/corncob/infiles/plasmid_ct_al_10_each.csv')
    plas_cnts.columns = plas_cnts.columns.to_series().str.replace('-','~')
    phag_cnts = read_csv('Data/corncob/infiles/phage_ct_al_10_each.csv')
    phag_cnts.columns = phag_cnts.columns.to_series().str.replace('-','~')


    #   PHAGES FIGURE
    phag_PC_CONT = read_csv('Data/corncob/results/PC_vs_CONT_phages.csv')
    phag_PC_CONT_age_gend = read_csv('Data/corncob/results/PC_vs_CONT_phages_cov_age_gender.csv')
    phag_PC_CONT_all_cov_imp = read_csv('Data/corncob/results/PC_vs_CONT_phages_all_covar_imp.csv')

    phag_PC_CONT_fdr = read_csv('Data/corncob/results/PC_vs_CONT_phages_fdrs.csv')
    phag_PC_CONT_age_gend_fdr = read_csv('Data/corncob/results/PC_vs_CONT_phages_fdrs_cov_age_gender.csv')
    phag_PC_CONT_all_cov_imp_fdr = read_csv('Data/corncob/results/PC_vs_CONT_phages_fdrs_all_covar_imp.csv')

    fdrs = [phag_PC_CONT_fdr, phag_PC_CONT_age_gend_fdr]
    dfs = [phag_PC_CONT, phag_PC_CONT_age_gend]
    res = []
    for df,name, fdr in zip(dfs, ['p_vs_c', 'p_vs_c_age_gend'],fdrs ):
        temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
        fdr.columns = ['otu', f'q_{name}']
        temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
        temp['stage_vs_0'] = name
        res.append(temp)

    name = 'p_vs_c_age_gend_all_cov_imp'
    df = phag_PC_CONT_all_cov_imp.copy()
    fdr = phag_PC_CONT_all_cov_imp_fdr.copy()
    temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
    fdr.columns = ['otu', f'q_{name}']
    temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
    temp['stage_vs_0'] = name
    res.append(temp)

    df1 = res[0].loc[res[0]['q_p_vs_c'] < 0.1][['OTU', 'stage_vs_0', 'Estimate', 'Pr(>|t|)']]

    cor_df =get_correl_with_host(df1, phag_cnts, mts_log)

    best_corr_per_mge = DataFrame()
    for mge in cor_df['mge'].unique():
        mge_df = cor_df.loc[cor_df['mge'] == mge].sort_values('q_val')
        best_corr_per_mge = concat([best_corr_per_mge,mge_df.iloc[0,:]], axis=1)
    best_corr_per_mge = best_corr_per_mge.T
    ##### Make Volcano plot

    col_dict = {'Bacteroides coprocola [ref_mOTU_v25_11279]': 'slateblue',
                'Bacteroides dorei/vulgatus [ref_mOTU_v25_02367]' : 'skyblue',
                'Bacteroides plebeius [ref_mOTU_v25_05069]' : 'darkblue',
                'Clostridium sp. CAG:7 [ref_mOTU_v25_07647]': 'lightpink',
                'Clostridiales species incertae sedis [meta_mOTU_v25_12452]' : 'firebrick',
                'Clostridiales species incertae sedis [meta_mOTU_v25_14513]' : 'firebrick',
                'Clostridiales species incertae sedis [meta_mOTU_v25_12473]' : 'firebrick',
                'Collinsella species incertae sedis [meta_mOTU_v25_12451]' : 'plum',
                'Enorma massiliensis [ref_mOTU_v25_02042]' : 'teal',
                'Fusobacterium nucleatum subsp. animalis [ref_mOTU_v25_01001]' : 'tomato',
                'Fusobacterium periodonticum [ref_mOTU_v25_00999]' : 'yellowgreen',
                'Prevotella copri [ref_mOTU_v25_03701]' : 'darkorange',
                'Prevotella species incertae sedis [meta_mOTU_v25_12765]' : (0.9921568627450981, 0.7490196078431373, 0.43529411764705883),
                'Streptococcus parasanguinis [ref_mOTU_v25_00312]' : (0.2, 0.6274509803921569, 0.17254901960784313),
                'Streptococcus sp. [ref_mOTU_v25_01350]' : 'gold',
                'Veillonella atypica [ref_mOTU_v25_01941]' : 'cornflowerblue',
                'Non Significant': 'gainsboro',}



    df5 = res[0].merge(best_corr_per_mge[['mge','bac']], left_on='OTU', right_on='mge', how='left')
    df5['-logp'] = -np.log10(df5['Pr(>|t|)'])
    df5.bac.fillna('Non Significant', inplace=True)
    plt.close('all')
    fig = plt.figure(figsize=(2.7,2.7))
    g = sns.scatterplot(x='Estimate', y='-logp', data = df5, hue=df5['bac'], \
                    palette = col_dict, linewidth=0, s=12, hue_order=col_dict.keys())
    plt.ylabel('')
    plt.xlabel('')

    plt.xlim(-7, 7)
    g.set(xticklabels=[])
    g.set(yticklabels=[])

    plt.hlines(3.78, -7,7, linestyle='dashed', colors='dimgrey', linewidth=1) # fdr 0.2

    g.get_legend().remove()

    plt.savefig(join(outdir, 'Fig3a_croncob.pdf'), dpi=300, bbox_inches='tight')


    ###PLASMID FIGURE

    plas_PC_CONT = read_csv('Data/corncob/results/PC_vs_CONT_plasmids.csv')
    plas_PC_CONT_age_gend = read_csv('Data/corncob/results/PC_vs_CONT_plasmids_cov_age_gender.csv')
    plas_PC_CONT_all_cov_imp = read_csv('Data/corncob/results/PC_vs_CONT_plasmids_all_covar_imp.csv')

    plas_PC_CONT_fdr = read_csv('Data/corncob/results/PC_vs_CONT_plasmids_fdrs.csv')
    plas_PC_CONT_age_gend_fdr = read_csv('Data/corncob/results/PC_vs_CONT_plasmids_fdrs_cov_age_gender.csv')
    plas_PC_CONT_all_cov_imp_fdr = read_csv('Data/corncob/results/PC_vs_CONT_plasmids_fdrs_all_covar_imp.csv')

    fdrs = [plas_PC_CONT_fdr, plas_PC_CONT_age_gend_fdr]
    dfs = [plas_PC_CONT, plas_PC_CONT_age_gend]
    res = []
    for df,name, fdr in zip(dfs, ['p_vs_c', 'p_vs_c_age_gend'],fdrs ):
        temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
        fdr.columns = ['otu', f'q_{name}']
        temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
        temp['stage_vs_0'] = name
        res.append(temp)

    name = 'p_vs_c_age_gend_all_cov_imp'
    df = plas_PC_CONT_all_cov_imp.copy()
    fdr = plas_PC_CONT_all_cov_imp_fdr.copy()
    temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
    fdr.columns = ['otu', f'q_{name}']
    temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
    temp['stage_vs_0'] = name
    res.append(temp)

    df1 = res[0].loc[res[0]['q_p_vs_c'] < 0.1][['OTU', 'stage_vs_0', 'Estimate', 'Pr(>|t|)']]

    cor_df =get_correl_with_host(df1, plas_cnts, mts_log)

    best_corr_per_mge = DataFrame()
    for mge in cor_df['mge'].unique():
        mge_df = cor_df.loc[cor_df['mge'] == mge].sort_values('q_val')
        best_corr_per_mge = concat([best_corr_per_mge,mge_df.iloc[0,:]], axis=1)
    best_corr_per_mge = best_corr_per_mge.T
    ##### Make Volcano plot


    df5 = res[0].merge(best_corr_per_mge[['mge','bac']], left_on='OTU', right_on='mge', how='left')
    df5['-logp'] = -np.log10(df5['Pr(>|t|)'])
    df5.bac.fillna('Non Significant', inplace=True)
    plt.close('all')
    fig = plt.figure(figsize=(2.7,2.7))
    g = sns.scatterplot(x='Estimate', y='-logp', data = df5, hue=df5['bac'], \
                    palette = col_dict, linewidth=0, s=12, hue_order=col_dict.keys())
    plt.ylabel('')
    plt.xlabel('')

    plt.xlim(-7, 7)
    g.set(xticklabels=[])
    g.set(yticklabels=[])

    plt.hlines(3.5, -7,7, linestyle='dashed', colors='dimgrey', linewidth=1) # fdr 0.2

    g.get_legend().remove()

    plt.savefig(join(outdir, 'Fig3b_croncob.pdf'), dpi=300, bbox_inches='tight')





if __name__ == "__main__":

    main()
