from pandas import read_pickle, concat,read_csv
from os.path import join
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
stats = importr('stats')
import warnings
warnings.filterwarnings("ignore")
from skbio.stats.distance import permanova


outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'

def main():
    g_meta = read_pickle('Data/metadata/full_german_pdac_metadata.df')
    s_meta = read_csv('Data/metadata/CNIO_metadata.tsv', sep='\t')
    g_meta.set_index('sample_alias', inplace=True)
    s_meta.set_index('sample_alias', inplace=True)
    s_meta = s_meta.reset_index().drop_duplicates(subset='sample_alias', keep='first').set_index('sample_alias')

    g_meta['study'] = 'German'
    s_meta['study'] = 'Spanish'

    meta = concat([g_meta[['subject_disease_status', 'study']],s_meta[['subject_disease_status', 'study']] ])
    meta = meta.loc[~(meta.subject_disease_status == 'NonElegible')]

    g_plas = read_pickle('Data/mOTUs/german_plasmid_data.pkl')
    g_phag = read_pickle('Data/mOTUs/german_phage_data.pkl')
    s_plas = read_pickle('Data/mOTUs/spanish_plasmid_data.pkl')
    s_phag = read_pickle('Data/mOTUs/spanish_phage_data.pkl')

    phage_df = g_phag.merge(s_phag, left_index=True, right_index=True, how='outer').fillna(0)
    plas_df = g_plas.merge(s_plas, left_index=True, right_index=True, how='outer').fillna(0)

    beta_div_metric = 'braycurtis'
    xa, ya = 'PC1', 'PC2'


    ordf = phage_df.T
    ph_mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat_ph = pcoa(ph_mat)
    df = pcoa_mat_ph.samples.copy()
    df.index = ordf.index
    x2ph = df.merge(meta, left_index=True, right_index=True, how = 'left')

    ordf = plas_df.T
    xa, ya = 'PC1', 'PC2'
    pl_mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat_pl = pcoa(pl_mat)
    df = pcoa_mat_pl.samples.copy()
    df.index = ordf.index
    x2pl = df.merge(meta, left_index=True, right_index=True, how = 'left')


    plt.close('all')

    fig, axes = plt.subplots(1, 2, sharex=True, figsize=(10,5))
    fig.suptitle('MGEs across German and Spanish studies')
    markers = {"German": "s", "Spanish": "X"}
    g = sns.scatterplot(x = xa, y=ya, data = x2ph, hue = 'subject_disease_status',s=15, ax=axes[0],
                        palette=dict(PC="salmon", CTR="darkcyan"), style='study', markers=markers)
    g.get_legend().remove()
    g = sns.scatterplot(x = xa, y=ya, data = x2pl, hue = 'subject_disease_status',s=15, ax=axes[1],
                        palette=dict(PC="salmon", CTR="darkcyan"), style='study', markers=markers)

    sns.move_legend(axes[1], "upper left", bbox_to_anchor=(1, 1))

    axes[0].set_xlabel(f'%s (%.2f%%)' % (xa, pcoa_mat_ph.proportion_explained.loc[xa]*100))
    axes[0].set_ylabel(f'%s (%.2f%%)' % (xa, pcoa_mat_ph.proportion_explained.loc[ya]*100))
    axes[1].set_xlabel(f'%s (%.2f%%)' % (xa, pcoa_mat_pl.proportion_explained.loc[xa]*100))
    axes[1].set_ylabel(f'%s (%.2f%%)' % (xa, pcoa_mat_pl.proportion_explained.loc[ya]*100))

    axes[0].set_title('Phages')
    axes[1].set_title('Plasmids')

    fig_name = 'FigS3_spa_germ_combined_PCAs.pdf'
    plt.savefig(join(outdir, fig_name), dpi=900, bbox_inches='tight')


    permanova_res = permanova(ph_mat, x2ph['study'], permutations=9999)
    print( 'phage study permanova pval = ' + str(permanova_res[5]))
    permanova_res = permanova(ph_mat, x2ph['subject_disease_status'], permutations=9999)
    print( 'phage disease pval = ' + str(permanova_res[5]))

    permanova_res = permanova(pl_mat, x2pl['study'], permutations=9999)
    print( 'plasmid study permanova pval = ' + str(permanova_res[5]))
    permanova_res = permanova(pl_mat, x2pl['subject_disease_status'], permutations=9999)
    print( 'plasmid disease pval = ' + str(permanova_res[5]))

if __name__ == "__main__":

    main()
