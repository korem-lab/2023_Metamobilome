import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from os.path import  join
from pandas import read_pickle, DataFrame, read_excel, read_csv
from mne.stats import fdr_correction
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

vsgv_file = 'Data/SGVFinder/vsgv_max_spacing_1_del_cooc_0.25msc_75.df'
meta_file = 'Data/metadata/CNIO_metadata.tsv'
icra_ab_file = 'Data/ICRA/ICRA_otu_ab_CNIO.df'
outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def get_sign_from_new_effect_size(num):
    if num > 0.5:
        yield 1
    if num < 0.5:
        yield -1

def assosciate_variable_region_and_cancer_mw(sample_region, meta, label_col, cancer_lbl, control_label, numbers_name , thresh):
    res = {}
    pc_list = meta.loc[meta[label_col] == cancer_lbl].index.to_list()
    cnt_list = meta.loc[meta[label_col] == control_label].index.to_list()

    for region in sample_region.columns:

        sub = sample_region.dropna(subset=[region])
        temp_pc = sub.loc[sub.index.isin(pc_list), region].to_list()
        temp_ctr = sub.loc[sub.index.isin(cnt_list), region].to_list()

        if len(temp_pc) < thresh or len(temp_ctr) <thresh:
            continue

        stat, pvalue = mannwhitneyu(temp_pc,temp_ctr, alternative='two-sided')
        effect_size = stat/(len(temp_pc)*len(temp_ctr))
        nums = f'({len(temp_pc)}|{len(temp_ctr)})'
        res[region] = [pvalue, effect_size, nums]

    res_df=DataFrame.from_dict(res, orient='index', columns=['pvalue', 'effect_size', numbers_name])
    res_df["qvalue"] = fdr_correction(res_df['pvalue'], alpha=0.1, method='indep')[1]
    return res_df

def main():
    #### CRC French cohort
    french_meta = read_pickle('Data/metadata/Frech_CRC_clean.pkl')
    c_vsgv_french = read_pickle('Data/SGVFinder/vsgv_bycnio_french_crc.df')

    c_vsgv_french.index = c_vsgv_french.index.to_series().str.replace('_united','')
    french_meta = french_meta.loc[french_meta.index.isin(c_vsgv_french.index)]
    french_norm_vs_cancer = assosciate_variable_region_and_cancer_mw(c_vsgv_french, french_meta,'Diagnosis', 'Cancer', 'Normal', 'numbers_F' , 20)
    mw_sig_french_crc = french_norm_vs_cancer.loc[french_norm_vs_cancer['pvalue'] < 0.05]

    #### CRC GUT cohort
    meta_crc_gut = read_excel('Data/metadata/CRC_Gut_metadata.xlsx')
    meta_gut = meta_crc_gut.loc[:,meta_crc_gut.nunique() > 1]
    meta_gut.drop(['Alias', 'AvgSpotLen', 'BioSample', 'Experiment', 'INSDC_center_name','INSDC_first_public',\
                  'INSDC_last_update','MBases','DATASTORE_provider','DATASTORE_region',\
                      'MBytes', 'SRA_accession',\
                      'DATASTORE_filetype','Library_Name', 'LoadDate'], axis=1 , inplace=True)
    meta_gut = meta_gut.set_index('Run')
    c_vsgv_gut = read_pickle('Data/SGVFinder/vsgv_bycnio_GUT_crc.df')
    crc_gut_mw_case_vs_control= assosciate_variable_region_and_cancer_mw(c_vsgv_gut, meta_gut, 'config','case', 'control', 'numbers_Gut' , 10)
    gut_sigs_wm = crc_gut_mw_case_vs_control.loc[crc_gut_mw_case_vs_control['pvalue'] < 0.05]

    ### CRC Nature Comm

    meta_NC = read_excel('Data/metadata/CRC_nature_com_metadata.xlsx')
    meta_NC=meta_NC.dropna(subset=['Alias'])
    meta_NC.title = meta_NC.title.replace('Stool sample from carcinoma', 'carcinoma')
    meta_NC.title = meta_NC.title.replace('Stool sample from advanced adenoma', 'adenoma')
    meta_NC.title = meta_NC.title.replace('Stool sample from controls', 'control')
    meta_NC.Alias = meta_NC.Alias.fillna(0).astype(int).astype(str)
    meta_NC = meta_NC.groupby(meta_NC.Alias).first()
    meta_NC.loc[:,meta_NC.nunique() > 1]

    c_vsgv_NC = read_pickle('Data/SGVFinder/vsgv_bycnio_natcom_crc.df')
    c_vsgv_NC.index = c_vsgv_NC.index.to_series().str.replace('_united','')
    NC_mw_control_vs_carcinoma = assosciate_variable_region_and_cancer_mw(c_vsgv_NC, meta_NC,'title', 'carcinoma', 'control', 'numbers_NC' , 10)
    NC_mw_sigs = NC_mw_control_vs_carcinoma.loc[NC_mw_control_vs_carcinoma['pvalue'] < 0.05]

    ###CRC PlosOne
    meta_po = read_pickle('Data/metadata/CRC_Plos_metadata.pkl')
    meta_po = meta_po.loc[:,meta_po.nunique() > 1]
    meta_po.drop(['/', 'AvgSpotLen', 'BioSample', 'Experiment', 'INSDC_last_update','MBases',\
                  'MBytes', 'sample_acc', 'Sample Name', 'sample_name', 'SRA_accession',\
                  'DATASTORE filetype','LibraryLayout'], axis=1, inplace=True)

    meta_po = meta_po.groupby(meta_po.Alias).first()

    meta_po.loc[meta_po.casectl == 1, "Diagnosis"] = 'cancer'
    meta_po.loc[meta_po.casectl == 0, 'Diagnosis'] = 'normal'
    meta_po.sex = meta_po.sex.replace(1, 'male')
    meta_po.sex = meta_po.sex.replace(2, 'female')
    meta_po.race = meta_po.race.replace(1, 'white')
    meta_po.race = meta_po.race.replace(2, 'black')
    meta_po.race = meta_po.race.replace(3, 'hispanic')
    meta_po.race = meta_po.race.replace(4, 'american indian/alaskan white')
    meta_po.race = meta_po.race.replace(5, 'asian or pacific islander')
    meta_po.race = meta_po.race.replace('.', np.nan)
    meta_po.race = meta_po.race.replace(8, np.nan)
    meta_po.race = meta_po.race.replace(9, np.nan)
    plos_vsgv = read_pickle('Data/SGVFinder/vsgv_bycnio_plosone_crc.df')
    plos_vsgv.index = plos_vsgv.index.to_series().str.replace('-27-0-0','')
    plos_one_wm_res = assosciate_variable_region_and_cancer_mw(plos_vsgv, meta_po, 'Diagnosis','cancer', 'normal', 'numbers_PO' , 10)
    plos_one_mw_sigs = plos_one_wm_res.loc[plos_one_wm_res['pvalue'] < 0.05]


    ### PDAC_German
    german_meta = read_pickle('Data/metadata/full_german_pdac_metadata.df')
    german_meta = german_meta.set_index('run_accession')
    german_vsgv_c = read_pickle('Data/SGVFinder/german_pdac_vsgv_bycnio.df')
    german_pdac_mw_res = assosciate_variable_region_and_cancer_mw(german_vsgv_c, german_meta,'subject_disease_status','PC', 'CTR', 'numbers_GP', 5)
    pdac_german_signif = german_pdac_mw_res.loc[german_pdac_mw_res['pvalue'] < 0.05]

    ###PDAC_Japan
    j_meta = read_pickle('Data/metadata/clean_japanese_pdac_metadata.df')
    j_meta = j_meta.set_index('Run')
    vsgv_j = read_pickle("Data/SGVFinder/japan_pdac_vsgv_bycnio.df")
    japan_pdac_mw_res = assosciate_variable_region_and_cancer_mw(vsgv_j, j_meta, 'Host_disease', 'PDAC', 'Control', 'numbers_JP', 10)
    japan_pdac_sigs = japan_pdac_mw_res.loc[japan_pdac_mw_res['pvalue'] < 0.05]

    ###CNIO

    meta = read_csv(meta_file, sep='\t')
    meta['source'] = meta.sample_alias.str.split('-').str[2]
    meta = meta.set_index('sample_alias')
    meta = meta[~meta.index.duplicated(keep='first')]
    meta_2 = meta.loc[meta['subject_disease_status'] != 'Pancreatitis']
    meta_2_stool = meta_2.loc[meta_2['source'] == 'ST']

    cnio_vsgv = read_pickle(vsgv_file)
    cnio_vsgv = cnio_vsgv.loc[cnio_vsgv.index.isin(meta_2_stool.index)]
    cnio_vsgv = cnio_vsgv.dropna(axis=1, how='all')
    cnio_res = assosciate_variable_region_and_cancer_mw(cnio_vsgv, meta_2_stool, 'subject_disease_status', 'PC', 'CTR', 'numbers_CNIO', 10)
    cnio_res = cnio_res.loc[cnio_res['qvalue'] < 0.2]

    ####combining
    cnio_res.rename(columns={'qvalue': 'cnio_q'}, inplace=True)
    cnio_res.loc[cnio_res['effect_size'] > 0.5, 'sign' ] = 1
    cnio_res.loc[cnio_res['effect_size'] < 0.5, 'sign' ] = -1
    mw_sig_french_crc.loc[mw_sig_french_crc['effect_size'] > 0.5, 'sign' ] = 1
    mw_sig_french_crc.loc[mw_sig_french_crc['effect_size'] < 0.5, 'sign' ] = -1
    gut_sigs_wm.loc[gut_sigs_wm['effect_size'] > 0.5, 'sign' ] = 1
    gut_sigs_wm.loc[gut_sigs_wm['effect_size'] < 0.5, 'sign' ] = -1
    plos_one_mw_sigs.loc[plos_one_mw_sigs['effect_size'] > 0.5, 'sign' ] = 1
    plos_one_mw_sigs.loc[plos_one_mw_sigs['effect_size'] < 0.5, 'sign' ] = -1
    NC_mw_sigs.loc[NC_mw_sigs['effect_size'] > 0.5, 'sign' ] = 1
    NC_mw_sigs.loc[NC_mw_sigs['effect_size'] < 0.5, 'sign' ] = -1
    pdac_german_signif.loc[pdac_german_signif['effect_size'] > 0.5, 'sign' ] = 1
    pdac_german_signif.loc[pdac_german_signif['effect_size'] < 0.5, 'sign' ] = -1
    japan_pdac_sigs.loc[japan_pdac_sigs['effect_size'] > 0.5, 'sign' ] = 1
    japan_pdac_sigs.loc[japan_pdac_sigs['effect_size'] < 0.5, 'sign' ] = -1
    cnio_res['CNIO_factor'] = -np.log10(cnio_res['pvalue'])* cnio_res['sign']
    mw_sig_french_crc['CRC_French_factor'] = -np.log10(mw_sig_french_crc['pvalue'])* mw_sig_french_crc['sign']
    pdac_german_signif['PDAC_Germany_factor'] = -np.log10(pdac_german_signif['pvalue'])* pdac_german_signif['sign']
    japan_pdac_sigs['PDAC_Japan_factor'] = -np.log10(japan_pdac_sigs['pvalue'])* japan_pdac_sigs['sign']
    merged = cnio_res.merge(mw_sig_french_crc, left_index=True, right_index=True, how='left')
    merged_crash = merged.merge(japan_pdac_sigs, left_index=True, right_index=True, how='left')
    merged_crash.drop(['pvalue_x', 'effect_size_x', 'sign_x', 'pvalue_y', 'effect_size_y'], axis=1, inplace=True)
    merged1 = merged_crash.merge(pdac_german_signif, left_index=True, right_index=True, how='left')
    gut_sigs_wm['CRC_GUT_factor'] = -np.log10(gut_sigs_wm['pvalue'])* gut_sigs_wm['sign']
    merged1.drop(['qvalue_x'], axis=1, inplace=True)
    merged2 = merged1.merge(gut_sigs_wm, left_index=True, right_index=True, how='left')
    NC_mw_sigs['CRC_NatCom_factor'] = -np.log10(NC_mw_sigs['pvalue'])* NC_mw_sigs['sign']
    merged2.drop(['sign_x', 'sign_x', 'pvalue_x','effect_size_x'], axis=1, inplace=True)
    merged3 = merged2.merge(NC_mw_sigs, left_index=True, right_index=True, how='left')
    merged3.drop(['sign_y', 'qvalue_y', 'pvalue_y', 'effect_size_y' , 'qvalue_x', 'sign_y','pvalue_x', 'effect_size_x', 'qvalue_y', 'sign_x'],axis=1, inplace=True)
    plos_one_mw_sigs['CRC_PO_factor'] = -np.log10(plos_one_mw_sigs['pvalue'])* plos_one_mw_sigs['sign']
    merged4 = merged3.merge(plos_one_mw_sigs, left_index=True, right_index=True, how='left')
    merged4['genome'] = merged4.index.to_series().str.split(':').str[0]
    merged4 = merged4.sort_values('genome')

    plt.close('all')
    fig = plt.figure(figsize = (6,2))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    ax= sns.heatmap(merged4[[ 'CRC_French_factor','CRC_GUT_factor' , 'CRC_PO_factor','CRC_NatCom_factor','PDAC_Germany_factor', 'PDAC_Japan_factor', 'CNIO_factor']].T, \
                yticklabels=False,  xticklabels=False,
                cmap = cmap, vmin = -5, vmax = 5)
    ax.axhline(y=0, color='k',linewidth=2) # works - top horizontal
    ax.axhline(y=7, color='k',linewidth=2)
    ax.axvline(x=0, color='k',linewidth=2)
    ax.axvline(x=29, color='k',linewidth=2)

    plt.savefig(join(outdir, 'fig1B_part1.pdf'),bbox_inches='tight', dpi=600)

    #### adding microbial abundance panel

    icra = read_pickle(icra_ab_file).T
    icra = icra.fillna(0)
    icra_meta = icra.merge(meta_2_stool, left_index=True, right_index=True, how='right')


    lst = ['445970.PRJNA19655', '1235788.PRJNA175977', '657318.PRJNA39159',
           '1121115.PRJNA195783', '748224.PRJNA46809', '411479.PRJNA18195']
    res =[]
    for genome in lst:

        res.append([genome, mannwhitneyu(icra_meta.loc[icra_meta['subject_disease_status'] == 'PC', genome].dropna(),
                                        icra_meta.loc[icra_meta['subject_disease_status'] == 'CTR', genome].dropna())[1]])

    res2 = DataFrame(res, columns=['genome','p'])
    res2['log_p'] = -np.log10(res2['p'])
    res2['factor'] = res2.log_p * -1

    fig = plt.figure(figsize = (4.8,0.26))
    df = merged4.merge(res2, left_on='genome', right_on='genome', how='left')
    ax= sns.heatmap(df[['factor']].T, \
                yticklabels=False,  xticklabels=False,
                cmap = cmap, vmin = -5, vmax = 5,
                    cbar=False)

    plt.savefig(join(outdir,'fig1B_part2.pdf'),bbox_inches='tight', dpi=600)



if __name__ == "__main__":
    main()
