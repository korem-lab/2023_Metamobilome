
from pandas import read_pickle, DataFrame, concat, read_csv
from os.path import join
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from scipy.stats import mannwhitneyu
from skbio.stats.ordination import pcoa
from scipy.spatial import distance
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
stats = importr('stats')
import warnings
warnings.filterwarnings("ignore")
from skbio.stats.distance import anosim, permanova

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'

def read_german_pdac_meta(metafile):
    meta = read_pickle(metafile)
    meta = meta.set_index('sample_alias')
    return meta

def norm_by_percent_alnd(pi_res_file, alnd_per_sample_file,thresh):
    pi_res = read_pickle(pi_res_file)
    pi_res = pi_res.truediv(pi_res.sum())

    alned_per_samp = read_pickle(alnd_per_sample_file)
    reads_stats = alned_per_samp.to_frame(name='alnd_rds')
    reads_stats['rds_counts'] = 1e7
    reads_stats['percent_alnd'] = reads_stats['alnd_rds']/reads_stats['rds_counts']

    for c in pi_res.columns:
        pi_res[c] = pi_res[c] * reads_stats.loc[c, 'percent_alnd']
    pi_res[pi_res < thresh] = 0
    pi_res = pi_res.fillna(0)
    return pi_res

def do_pca(data,meta, meta_vars, figure_out_path, figure_name_string):

    beta_div_metric = 'braycurtis'
    ordf = data.T
    xa, ya = 'PC1', 'PC2'
    mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat = pcoa(mat)
    df = pcoa_mat.samples.copy()
    df.index = ordf.index
    x2 = df.merge(meta, left_index=True, right_index=True, how = 'left')
    print(figure_name_string)
    print(xa, pcoa_mat.proportion_explained.loc[xa]*100)
    print(ya, pcoa_mat.proportion_explained.loc[ya]*100)

    for variable in meta_vars:
        plt.close('all')
        plt.figure(figsize=(1.8,1.8))
        g = sns.scatterplot(x = xa, y=ya, data = x2.round(2), hue = variable, palette=dict(PC="salmon", CTR="darkcyan").values(),s=12)
        g.get_legend().remove()

        plt.xlabel("")
        plt.ylabel("")
        g.set(xticklabels=[])
        g.set(yticklabels=[])

        fig_name = 'Fig2_'+figure_name_string + '_PCA_by_' + beta_div_metric + '_' + variable + '.pdf'
        plt.savefig(join(figure_out_path, fig_name), dpi=300)
        permanova_res = permanova(mat, x2[variable], permutations=9999)
        print( 'permanova pval = ' + str(permanova_res[5]))


    print(figure_name_string)
    print(xa, pcoa_mat.proportion_explained.loc[xa]*100)
    print(ya, pcoa_mat.proportion_explained.loc[ya]*100)


def find_distnace_from_centroid(df):
    dists = []
    centroid = df.mean().to_numpy()
    for i in df.index:
        dists.append(distance.braycurtis(centroid, df.loc[i,:].values))
    return dists

def distance_from_centroid(data, meta, figure_out_path, figure_name_string):
    data_meta = data.T.merge(meta[['subject_disease_status']], left_index=True, right_index=True, how='left')
    sub_pc = data_meta.loc[data_meta.subject_disease_status == 'PC']
    sub_ctr = data_meta.loc[data_meta.subject_disease_status == 'CTR']
    sub_pc.drop(columns=['subject_disease_status'], inplace=True)
    sub_ctr.drop(columns=['subject_disease_status'], inplace=True)
    pc_dists = find_distnace_from_centroid(sub_pc)
    ctr_dists = find_distnace_from_centroid(sub_ctr)

    temp = concat([DataFrame([pc_dists, ['PC']*len(pc_dists)]).T,DataFrame([ctr_dists, ['CTR']*len(ctr_dists)]).T])

    plt.close('all')
    plt.figure(figsize=(0.8,1.8))
    g = sns.boxplot(x = 1,y = 0, data = temp, palette=dict(PC="salmon", CTR="darkcyan"),showfliers = False)
    sns.swarmplot(x = 1,y = 0, data = temp, color='black', s=1.8)

    plt.xlabel("")
    plt.ylabel("")
    plt.ylim(0.4,1.2)
    g.set(xticklabels=[])
    fig_name = 'Fig2_' + figure_name_string + 'german_dist_from_centroid.pdf'
    plt.savefig(join(figure_out_path, fig_name), dpi=300)
    print(figure_name_string + 'dist_from_centroid.pdf')
    print(mannwhitneyu(temp.loc[temp[1] == 'PC', 0].to_list(),temp.loc[temp[1] == 'CTR', 0].to_list(), alternative='two-sided')[1])


def number_of_plasmids_per_sample(pres_abs, meta, figure_out_path, figure_name_string):
    plas_per_sample = pres_abs.sum(0).rename('no_of_diff_elements')
    plas_sum_df = meta.merge(plas_per_sample,right_index=True, left_index=True,  how='right' )
    plas_sum_df = plas_sum_df.sort_values('subject_disease_status', ascending=False)
    plt.close('all')
    plt.figure(figsize=(0.8,1.8))
    ax = plt.subplot(111)

    g = sns.boxplot(x = 'subject_disease_status', y = 'no_of_diff_elements', ax = ax, data = plas_sum_df, palette=dict(PC="salmon", CTR="darkcyan"),showfliers = False)
    sns.swarmplot(x = 'subject_disease_status',y = 'no_of_diff_elements', data = plas_sum_df, color='black', s=1.8)

    plt.xlabel("")
    plt.ylabel("")
    plt.ylim(0,15000)
    g.set(xticklabels=[])
    fig_name = 'Fig2_no_of_' + figure_name_string + '_per_sample_by_disease_state.pdf'
    plt.savefig(join(figure_out_path, fig_name), dpi=300)
    print(figure_out_path,figure_name_string)
    print(mannwhitneyu(plas_sum_df.loc[plas_sum_df['subject_disease_status'] == 'PC', 'no_of_diff_elements'].values,
                                         plas_sum_df.loc[plas_sum_df['subject_disease_status'] == 'CTR', 'no_of_diff_elements'].values, alternative='two-sided')[1])



def run_analysis(pi_res, meta, meta_vars, figures_path, figure_name_string , pres_abs):

    do_pca(pi_res,meta, meta_vars, figures_path, figure_name_string)
    distance_from_centroid(pi_res, meta, figures_path, figure_name_string)
    number_of_plasmids_per_sample(pres_abs, meta, figures_path, figure_name_string)




def main(mge_info_file, meta_file, pi_res_file, alnd_per_sample_file, pi_thresh,figures_path):
    mmge_info = read_csv(mge_info_file, sep=',')
    meta_vars = ['subject_disease_status']

    meta = read_german_pdac_meta(meta_file)
    pi_res = norm_by_percent_alnd(pi_res_file, alnd_per_sample_file, pi_thresh)
    merged = pi_res.T.merge(meta, left_index=True, right_on='run_accession', how='left')
    pi_res = merged.loc[merged.subject_disease_status != 'Pancreatitis']
    pi_res = pi_res.drop(meta.columns, axis=1)
    pi_res = pi_res.T
    pi_res = pi_res.loc[(pi_res != 0).sum(1) > 5]
    pi_res = pi_res.merge(mmge_info[['MGEs_id', 'status']], left_index=True, right_on='MGEs_id', how='left').set_index(['MGEs_id','status'])
    pi_res = pi_res.swaplevel('MGEs_id', 1)


    pres_abs = (pi_res > 0) * 1
    pres_abs = pres_abs[pres_abs.sum(1) > 0]
    run_analysis(pi_res.loc['phage'], meta, meta_vars, figures_path, 'German_Stool_Phages_samples', pres_abs.loc['phage'])
    run_analysis(pi_res.loc['plasmid'], meta, meta_vars, figures_path, 'German_Stool_Plasmid_samples', pres_abs.loc['plasmid'])

if __name__ == "__main__":

    main('Data/random/mMGE_metadata.csv',
         'Data/metadata/full_german_pdac_metadata.df',
         'Data/ICRA/german_pdac_mMGE100k_cover_res_ra.pkl',
         'Data/ICRA/german_pdac_mMGE100k_alnd_rds_per_sample.pkl',
         0.0000001,outdir)







