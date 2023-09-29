from statistics import median, mean
import pandas as pd
from pandas import read_pickle, DataFrame, to_pickle, concat, read_csv, isnull, Series, read_excel
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity

import warnings
warnings.filterwarnings("ignore")
from skbio.diversity import alpha

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def add_p_val(lft,rgt,y,h,p):
    plt.plot([lft, lft, rgt, rgt], [y, y+h, y+h, y], lw=1.5, c='k')
    plt.text((lft + rgt) * .5, y+h, ('n.s.' if p > 0.15 else 'p < %.2g' if p > 0.001 else 'p < %.1g') % max(p+1e-20, 1e-20), ha='center', va='bottom', color='k')


def read_meta_cnio(meta_file):
    meta = read_csv(meta_file, sep='\t')
    meta['source'] = meta.sample_alias.str.split('-').str[2]
    meta = meta.set_index('sample_alias')
    meta = meta[~meta.index.duplicated(keep='first')]
    meta = meta.loc[meta.subject_disease_status != 'NonElegible']
    meta = meta.dropna(subset=['diabin'])
    meta.diabin = meta.diabin.astype('int')
    return meta

def make_fig(shannon_df, col_dict, fig_type, cohort):
    plt.close('all')
    plt.figure(figsize=(2.5, 2))
    temp_df = shannon_df.loc[shannon_df.data_type == fig_type]
    temp_df.stage2 = temp_df.stage2.astype('str')
    g = sns.boxplot(x='stage2', y='Shannon Div', data=temp_df, palette=col_dict, showfliers = False, linewidth=0.75)
    sns.swarmplot(x='stage2', y='Shannon Div', data=temp_df, color='black', s=1.8)
    plt.ylim(1,14)
    plt.xlabel("")
    plt.ylabel("")
    g.set(xticklabels=[])
    g.set(yticklabels=[])
    figname = f'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/figS5_{cohort}_{fig_type}.pdf'
    plt.savefig(figname,dpi=300)


def do_pca_stage12_vs_cont(data,meta,figure_out_path, figure_name_string, fig_type):
    col_dict = {'0':"darkcyan", '1-2':"salmon"}
    beta_div_metric = 'braycurtis'
    ordf = data.T
    xa, ya = 'PC1', 'PC2'

    mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat = pcoa(mat)
    df = pcoa_mat.samples.copy()
    df.index = ordf.index
    x2 = df.merge(meta, left_index=True, right_index=True, how = 'left')
    x2['stage2'] = x2['stage2'].astype('str')
    bdiv_df = mat.to_data_frame()
    bdiv_df.columns = ordf.index
    bdiv_df.index = ordf.index

    print(fig_type)
    print(xa, pcoa_mat.proportion_explained.loc[xa]*100)
    print(ya, pcoa_mat.proportion_explained.loc[ya]*100)
    plt.close('all')
    plt.figure(figsize=(1.8,1.8))
    g = sns.scatterplot(x = xa, y=ya, data = x2.round(2), hue = 'stage2',palette = col_dict, s=12, hue_order=['0', '1-2']) #
    plt.xlabel("")
    plt.ylabel("")
    g.set(xticklabels=[])
    g.set(yticklabels=[])
    g.get_legend().remove()
    permanova_res_all = permanova(mat, x2['stage2'], permutations=9999)
    print( 'permanova all pval = ' + str(permanova_res_all[5]))
    figname = f'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/fig4_early_vs_cont_beta_div_{fig_type}.pdf'
    plt.savefig(figname,dpi=300)



def do_pca_stage(data,meta,figure_out_path, figure_name_string, fig_type):
    col_dict = {'0':"darkcyan", '1-2':"salmon", '3':"palevioletred", '4':"purple"}

    beta_div_metric = 'braycurtis'
    ordf = data.T
    xa, ya = 'PC1', 'PC2'
    mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat = pcoa(mat)
    df = pcoa_mat.samples.copy()
    df.index = ordf.index
    x2 = df.merge(meta, left_index=True, right_index=True, how = 'left')
    x2['stage2'] = x2['stage2'].astype('str')
    bdiv_df = mat.to_data_frame()
    bdiv_df.columns = ordf.index
    bdiv_df.index = ordf.index
    print(fig_type)
    print(xa, pcoa_mat.proportion_explained.loc[xa]*100)
    print(ya, pcoa_mat.proportion_explained.loc[ya]*100)
    plt.close('all')
    plt.figure(figsize=(2,2))
    g = sns.scatterplot(x = xa, y=ya, data = x2.round(2), hue = 'stage2',palette = col_dict, s=35, hue_order=['0', '1-2', '3', '4']) #
    g.get_legend().remove()
    g.set(xlabel=None)
    g.set(ylabel=None)

    plt.xlabel("")
    plt.ylabel("")
    g.set(xticklabels=[])
    g.set(yticklabels=[])


    permanova_res_all = permanova(mat, x2['stage2'], permutations=9999)
    print( 'permanova all pval = ' + str(permanova_res_all[5]))
    figname = f'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/figS6_all_stage_pcoa_{fig_type}.pdf'
    plt.savefig(figname,dpi=300)
    flierprops = dict(markersize=2, linestyle='none')
    res = x2[['PC1', 'PC2', 'stage2']]
    plt.close('all')
    plt.figure(figsize=(2,0.4))
    g2 = sns.boxplot(y='stage2', x='PC1', data=res, palette = col_dict, order=['0', '1-2', '3', '4'], linewidth=0.5,flierprops=flierprops)
    g2 = sns.boxplot(y='stage2', x='PC1', data=res, palette = col_dict, order=['0','1-2', '3', '4'], linewidth=0.5,flierprops=flierprops)
    figname = f'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/figS6_beta_div_PC1_{fig_type}.pdf'
    g2.set(xticklabels=[])
    g2.set(yticklabels=[])
    g2.set(xlabel=None)
    g2.set(ylabel=None)
    g2.tick_params(left=False, bottom=False)
    plt.savefig(figname,dpi=300)

    plt.close('all')
    plt.figure(figsize=(0.4,2))
    g2 = sns.boxplot(x='stage2', y='PC2', data=res, palette = col_dict, order=['0','1-2', '3', '4'], linewidth=0.5,flierprops=flierprops)
    figname = f'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/figS6_beta_div_PC2_{fig_type}.pdf'
    g2.set(xticklabels=[])
    g2.set(yticklabels=[])
    g2.set(xlabel=None)
    g2.set(ylabel=None)
    g2.tick_params(left=False, bottom=False)
    plt.savefig(figname,dpi=300)


def main(outdir):
    # MAKE_FIG_WITH_LEGEND = True

    #importing CNIO DATA
    meta = read_meta_cnio('Data/metadata/CNIO_metadata.tsv')
    stool_samples = meta.loc[meta['source'] == 'ST'].index
    meta_stool = meta.loc[stool_samples]
    meta_stool.loc[meta_stool['stage'] == 1, 'stage2'] = '1-2'
    meta_stool.loc[meta_stool['stage'] == 2, 'stage2'] = '1-2'
    meta_stool.loc[meta_stool['stage'] == 3, 'stage2'] = '3'
    meta_stool.loc[meta_stool['stage'] == 4, 'stage2'] = '4'
    meta_stool['stage2'] = meta_stool['stage2'].fillna('0')

    plasmid_ra = read_pickle('Data/ICRA/CNIO_plasmid_ra_processed_100K_and_below.pkl')
    phage_ra = read_pickle('Data/ICRA/CNIO_phage_ra_processed_100K_and_below.pkl')
    motus = read_csv('Data/mOTUs/CNIO_motus2.tsv', sep='\t')
    motus = motus.loc[:, motus.columns != 'MMPC-3109233-ST-0']
    motus = motus.set_index('#consensus_taxonomy')
    motus = motus.fillna(0)
    motus.drop('-1', inplace=True)
    motus = motus.loc[:,motus.columns.isin(phage_ra.columns)]

    #collecting CNIO SHannon div
    shannon_df = DataFrame()
    for cur_df, label in zip([motus,plasmid_ra,phage_ra],['mts', 'plasmid','phage']):
        sdis = []
        for col in cur_df.columns:
            sdis.append(alpha.shannon(cur_df[col]))
        r = Series(sdis, index=cur_df.columns)
        r.name = 'Shannon Div'
        df_4fig = meta_stool[['subject_disease_status']].merge(DataFrame(r),right_index=True, left_index=True,  how='right' )

        df_4fig['data_type'] = label
        shannon_df = concat([shannon_df, df_4fig])
    shannon_df = shannon_df.loc[shannon_df['subject_disease_status'] != 'Pancreatitis']
    shannon_df = shannon_df.merge(meta_stool[['stage2']], left_index=True, right_index=True, how='left')
    shannon_df.loc[shannon_df['subject_disease_status'] == 'CTR', 'stage'] = 0
    shannon_df = shannon_df.sort_values('stage')

    col_dict = {'0':"darkcyan", '1-2':"salmon", '3':"palevioletred", '4':"purple"}

    for type in ['mts', 'plasmid','phage']:
        make_fig(shannon_df, col_dict, type, 'CNIO')


    #import German Data
    g_meta = pd.read_pickle('Data/metadata/full_german_pdac_metadata.df')
    g_meta = g_meta.set_index('sample_alias')
    g_meta.loc[g_meta['stage'] == 1, 'stage2'] = '1-2'
    g_meta.loc[g_meta['stage'] == 2, 'stage2'] = '1-2'
    g_meta.loc[g_meta['stage'] == 3, 'stage2'] = '3'
    g_meta.loc[g_meta['stage'] == 4, 'stage2'] = '4'
    g_meta['stage2'] = g_meta['stage2'].fillna('0')


    g_plasmid_ra = read_pickle('Data/ICRA/german_pdac_plasmid_ra_processed_100K_and_below.pkl')
    g_phage_ra = read_pickle('Data/ICRA/german_pdac_phage_ra_processed_100K_and_below.pkl')

    g_motus = read_csv('Data/mOTUs/German_PDAC_motus2.tsv', sep='\t')
    g_motus = g_motus.set_index('#consensus_taxonomy')
    g_motus = g_motus.fillna(0)
    g_motus.drop('-1', inplace=True)
    g_motus.columns = [x.split('_')[0] for x in g_motus.columns]
    g_motus = g_motus.T.merge(g_meta[['run_accession', 'alias' ]], left_index=True, right_on='run_accession', how='left').set_index('alias').drop('run_accession', axis=1).T


    shannon_df = DataFrame()
    for cur_df, label in zip([g_motus,g_plasmid_ra,g_phage_ra],['mts', 'plasmid','phage']):
        sdis = []
        for col in cur_df.columns:
            sdis.append(alpha.shannon(cur_df[col]))
        r = Series(sdis, index=cur_df.columns)
        r.name = 'Shannon Div'
        df_4fig = g_meta[['subject_disease_status']].merge(DataFrame(r),right_index=True, left_index=True,  how='right' )
        df_4fig['data_type'] = label
        shannon_df = concat([shannon_df, df_4fig])
    shannon_df = shannon_df.loc[shannon_df['subject_disease_status'] != 'Pancreatitis']
    shannon_df = shannon_df.merge(g_meta[['stage2']], left_index=True, right_index=True, how='left')
    shannon_df.loc[shannon_df['subject_disease_status'] == 'CTR', 'stage'] = 0
    shannon_df = shannon_df.sort_values('stage2')
    shannon_df.loc[shannon_df['subject_disease_status'] == 'CTR', 'stage'] = 0
    shannon_df.dropna(subset=['stage2'], inplace=True)

    for type in ['mts', 'plasmid','phage']:
        make_fig(shannon_df, col_dict, type, 'German')

    for df, fig_type in zip([motus, plasmid_ra, phage_ra],['s_mts', 's_plasmid','s_phage']):
        do_pca_stage(df, meta_stool, 'bla', 'bla', fig_type)

    for df, fig_type in zip([g_motus,g_plasmid_ra,g_phage_ra],['g_mts', 'g_plasmid','g_phage']):
        do_pca_stage(df, g_meta, 'bla', 'bla', fig_type)

    for df, fig_type in zip([motus, plasmid_ra, phage_ra],['s_mts', 's_plasmid','s_phage']):
        do_pca_stage12_vs_cont(df, meta_stool, 'bla', 'bla', fig_type)

    for df, fig_type in zip([g_motus,g_plasmid_ra,g_phage_ra],['g_mts', 'g_plasmid','g_phage']):
        do_pca_stage12_vs_cont(df, g_meta, 'bla', 'bla', fig_type)


if __name__ == "__main__":
    main(outdir)
