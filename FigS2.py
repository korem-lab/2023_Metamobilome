from pandas import read_csv,DataFrame,concat,read_excel, read_pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from skbio.stats.ordination import pcoa

from skbio.stats.distance import permanova
from skbio.diversity import beta_diversity
from os.path import join

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def main():

    ########pcoa of bacterial species across cohorts
    #### CRC French
    french_meta = read_pickle('Data/metadata/Frech_CRC_clean.pkl')
    french_otu = read_pickle('Data/mOTUs/french_mOTUs2.df')
    french_otu.columns = french_otu.columns.droplevel([1,2])
    french_otu = french_otu.drop('-1', axis=1)
    french_otu = french_otu.truediv(french_otu.sum(1), axis=0)
    french_otu.index = french_otu.index.to_series().str.replace('_united','')
    french_meta = french_meta.loc[french_meta.index.isin(french_otu.index)]
    french_df = french_otu.merge(french_meta, left_index=True, right_index=True, how='left')
    french_df = french_df.loc[french_df.Diagnosis.isin(['Cancer', 'Normal'])]
    french_df['phenotype'] = np.where(french_df['Diagnosis'] == 'Normal', 'Control', 'Cancer')


    #### CRC GUT cohort
    meta_crc_gut = read_excel('Data/metadata/CRC_Gut_metadata.xlsx')
    gut_meta = meta_crc_gut.loc[:,meta_crc_gut.nunique() > 1]
    gut_meta.drop(['Alias', 'AvgSpotLen', 'BioSample', 'Experiment', 'INSDC_center_name','INSDC_first_public',\
                  'INSDC_last_update','MBases','DATASTORE_provider','DATASTORE_region',\
                      'MBytes', 'SRA_accession',\
                      'DATASTORE_filetype','Library_Name', 'LoadDate'], axis=1 , inplace=True)
    gut_meta = gut_meta.set_index('Run')
    gut_otu = read_pickle('Data/mOTUs/GUT_mOTUs2.df')
    gut_otu.columns = gut_otu.columns.droplevel([1,2])
    gut_otu = gut_otu.drop('-1', axis=1)
    gut_otu = gut_otu.truediv(gut_otu.sum(1), axis=0)
    gut_df = gut_otu.merge(gut_meta, left_index=True, right_index=True, how='left')
    gut_df['phenotype'] = np.where(gut_df['config'] == 'control', 'Control', 'Cancer')

    ### CRC Nature Comm
    NC_meta = read_excel('Data/metadata/CRC_nature_com_metadata.xlsx')
    NC_meta=NC_meta.dropna(subset=['Alias'])
    NC_meta.title = NC_meta.title.replace('Stool sample from carcinoma', 'carcinoma')
    NC_meta.title = NC_meta.title.replace('Stool sample from advanced adenoma', 'adenoma')
    NC_meta.title = NC_meta.title.replace('Stool sample from controls', 'control')
    NC_meta.Alias = NC_meta.Alias.fillna(0).astype(int).astype(str)
    NC_meta = NC_meta.groupby(NC_meta.Alias).first()
    NC_meta.loc[:,NC_meta.nunique() > 1]
    NC_otu = read_pickle('Data/mOTUs/NC_mOTUs2.df')
    NC_otu.index = NC_otu.index.to_series().str.replace('_united','')
    NC_otu.columns = NC_otu.columns.droplevel([1,2])
    NC_otu = NC_otu.drop('-1', axis=1)
    NC_otu = NC_otu.truediv(NC_otu.sum(1), axis=0)

    NC_df = NC_otu.merge(NC_meta, left_index=True, right_index=True, how='left')
    NC_df = NC_df.loc[NC_df.title.isin(['carcinoma', 'control'])]
    NC_df['phenotype'] = np.where(NC_df['title'] == 'control', 'Control', 'Cancer')

    ###CRC PlosOne
    PO_meta = read_pickle('Data/metadata/CRC_Plos_metadata.pkl')
    PO_meta = PO_meta.loc[:,PO_meta.nunique() > 1]
    PO_meta.drop(['/', 'AvgSpotLen', 'BioSample', 'Experiment', 'INSDC_last_update','MBases',\
                  'MBytes', 'sample_acc', 'Sample Name', 'sample_name', 'SRA_accession',\
                  'DATASTORE filetype','LibraryLayout'], axis=1, inplace=True)


    PO_meta = PO_meta.groupby(PO_meta.Alias).first()

    PO_meta.loc[PO_meta.casectl == 1, "Diagnosis"] = 'cancer'
    PO_meta.loc[PO_meta.casectl == 0, 'Diagnosis'] = 'normal'
    PO_meta.sex = PO_meta.sex.replace(1, 'male')
    PO_meta.sex = PO_meta.sex.replace(2, 'female')
    PO_meta.race = PO_meta.race.replace(1, 'white')
    PO_meta.race = PO_meta.race.replace(2, 'black')
    PO_meta.race = PO_meta.race.replace(3, 'hispanic')
    PO_meta.race = PO_meta.race.replace(4, 'american indian/alaskan white')
    PO_meta.race = PO_meta.race.replace(5, 'asian or pacific islander')
    PO_meta.race = PO_meta.race.replace('.', np.nan)
    PO_meta.race = PO_meta.race.replace(8, np.nan)
    PO_meta.race = PO_meta.race.replace(9, np.nan)
    PO_otu = read_pickle('Data/mOTUs/PO_mOTUs2.df')
    PO_otu.index = PO_otu.index.to_series().str.replace('-27-0-0','')
    PO_otu.columns = PO_otu.columns.droplevel([1,2])
    PO_otu = PO_otu.drop('-1', axis=1)
    PO_otu = PO_otu.truediv(PO_otu.sum(1), axis=0)
    PO_otu = PO_otu.loc[PO_meta.dropna(subset=['Diagnosis']).index]
    PO_df = PO_otu.merge(PO_meta, left_index=True, right_index=True, how='left')
    PO_df['phenotype'] = np.where(PO_df['Diagnosis'] == 'normal', 'Control', 'Cancer')

    ### PDAC_German
    german_meta = read_pickle('Data/metadata/full_german_pdac_metadata.df')
    german_meta = german_meta.set_index('run_accession')
    german_otu = read_csv('Data/mOTUs/German_PDAC_motus2.tsv', sep='\t')#, skiprows=2)
    german_otu = german_otu.set_index('#consensus_taxonomy').T.iloc[:,:-1]
    german_otu = german_otu.truediv(german_otu.sum(1), axis=0)
    german_otu.columns = german_otu.columns.str.split(' ').str[-1].str.replace(']','').str.replace('[','')

    german_otu.index = german_otu.index.to_series().str.split('_').str[0]
    german_df = german_otu.merge(german_meta, left_index=True, right_index=True, how='left')
    german_df['phenotype'] = np.where(german_df['subject_disease_status'] == 'CTR', 'Control', 'Cancer')

    ###PDAC_Japan
    jap_meta = read_pickle('Data/metadata/clean_japanese_pdac_metadata.df')
    jap_meta = jap_meta.set_index('Run')
    jap_otu = read_csv("Data/mOTUs/Japan_PDAC_motus2.tsv", sep='\t')
    jap_otu = jap_otu.set_index('#consensus_taxonomy').T.iloc[:,:-1]
    jap_otu = jap_otu.truediv(jap_otu.sum(1), axis=0)
    jap_otu.columns = jap_otu.columns.str.split(' ').str[-1].str.replace(']','').str.replace('[','')

    jap_otu.index = jap_otu.index.to_series().str.split('_').str[0]
    jap_df = jap_otu.merge(jap_meta, left_index=True, right_index=True, how='inner')
    jap_df['phenotype'] = np.where(jap_df['Host_disease'] == 'Control', 'Control', 'Cancer')


    ###CNIO
    def read_meta_cnio(meta_file):
        meta = read_csv(meta_file, sep='\t')
        meta['source'] = meta.sample_alias.str.split('-').str[2]
        meta = meta.set_index('sample_alias')
        meta = meta[~meta.index.duplicated(keep='first')]
        meta = meta.loc[meta.subject_disease_status != 'NonElegible']
        meta = meta.dropna(subset=['diabin'])
        meta.intervention_full = meta.intervention_full.fillna('NONE')
        meta.diabin = meta.diabin.astype('int')
        return meta
    cnio_meta = read_meta_cnio('Data/metadata/CNIO_metadata.tsv')
    cnio_otu = read_csv('Data/mOTUs/CNIO_motus2.tsv',sep='\t')
    cnio_otu = cnio_otu.set_index('#consensus_taxonomy').T.iloc[:,:-1]
    cnio_otu = cnio_otu.truediv(cnio_otu.sum(1), axis=0)
    cnio_otu.columns = cnio_otu.columns.str.split(' ').str[-1].str.replace(']','').str.replace('[','')
    cnio_df = cnio_otu.merge(cnio_meta, left_index=True, right_index=True, how='left')
    cnio_df = cnio_df.loc[cnio_df['subject_disease_status'].isin(['PC','CTR'])].loc[cnio_df['source'] == 'ST']
    cnio_df['phenotype'] = np.where(cnio_df['subject_disease_status'] == 'CTR', 'Control', 'Cancer')


    dfs = [cnio_df, german_df, jap_df, PO_df, NC_df, french_df, gut_df]

    labels = ['PDAC Spain', 'PDAC Germany', 'PDAC Japan', 'CRC USA', 'CRC Austria', 'CRC France', 'CRC China']

    select_df = DataFrame()
    for df, lbl in zip(dfs, labels):
        df['cohort'] = lbl
        tmp = df[['cohort', 'phenotype']]
        select_df = concat([select_df, tmp])



    df1 = french_otu.T.merge(gut_otu.T, left_index=True, right_index=True, how='outer').merge(NC_otu.T, left_index=True, right_index=True, how='outer').merge(PO_otu.T, left_index=True, right_index=True, how='outer')

    merged_otu = df1.merge(german_otu.T, left_index=True, right_index=True, how='outer').merge(jap_otu.T, left_index=True, right_index=True, how='outer').merge(cnio_otu.T, left_index=True, right_index=True, how='outer')
    merged_otu = merged_otu.T.fillna(0).merge(select_df, left_index=True, right_index=True, how='inner')
    merge_otu_controls = merged_otu.loc[merged_otu['phenotype'] == 'Control'].drop(['phenotype', 'cohort'], axis=1)#,'445970.PRJNA19655', '1235788.PRJNA175977', '657318.PRJNA39159',
     #      '1121115.PRJNA195783', '748224.PRJNA46809', '411479.PRJNA18195'],axis=1)
    merge_otu_controls = merge_otu_controls.loc[:,(merge_otu_controls > 0).sum() > 1]

    my_palette = {"PDAC Spain":"crimson",
               "PDAC Germany":"pink",
               "PDAC Japan":"mediumblue",
               "CRC USA":"seagreen",
               "CRC Austria":"tomato",
               "CRC France":"coral",
               "CRC China":"cornflowerblue"}

    beta_div_metric = 'jaccard'#'jaccard'
    ordf = merge_otu_controls
    xa, ya = 'PC1', 'PC2'
    mat = beta_diversity(beta_div_metric, ordf.values)
    pcoa_mat = pcoa(mat)
    df = pcoa_mat.samples.copy()
    df.index = ordf.index
    x2 = df.merge(select_df, left_index=True, right_index=True, how = 'left')

    bdiv_df = mat.to_data_frame()
    bdiv_df.columns = ordf.index
    bdiv_df.index = ordf.index

    plt.close('all')
    plt.figure(figsize=(4,3))
    g = sns.scatterplot(x=xa,y=ya, data = x2, hue = 'cohort', s=14, alpha=0.6, palette=my_palette, linewidth=0,
                        hue_order = ["PDAC Spain", "PDAC Germany", "CRC Austria","CRC France","CRC USA","CRC China","PDAC Japan"])#, style="phenotype")
    plt.xlabel('%s (%.2f%%)' % (xa, pcoa_mat.proportion_explained.loc[xa]*100))
    plt.ylabel('%s (%.2f%%)' % (ya, pcoa_mat.proportion_explained.loc[ya]*100))
    sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))

    geo_dict = {'CRC France' : 'Europe',  'CRC Austria' : 'Europe', 'PDAC Germany': 'Europe', 'CRC China':'Asia', 'PDAC Japan' : 'Asia', 'CRC USA' : 'USA', 'PDAC Spain' : 'Europe' }
    x2['geo'] = x2['cohort'].map(geo_dict)
    permanova_res = permanova(mat, x2['geo'], permutations=9999)
    print( 'permanova pval = ' + str(permanova_res[5]))
    plt.savefig(join(outdir,'FigS2.pdf'), dpi=600, bbox_inches='tight')





if __name__ == "__main__":
    main()
