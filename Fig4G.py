from pandas import read_csv
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def main():


    #####MOTUS
    s0_2 = read_csv('Data/corncob/stage2_vs_CONT_motus.csv')
    s0_3 = read_csv('Data/corncob/stage3_vs_CONT_motus.csv')
    s0_4 = read_csv('Data/corncob/stage4_vs_CONT_motus.csv')

    s0_2_fdr = read_csv('Data/corncob/stage2_vs_CONT_motus_fdrs.csv')
    s0_3_fdr = read_csv('Data/corncob/stage3_vs_CONT_motus_fdrs.csv')
    s0_4_fdr = read_csv('Data/corncob/stage4_vs_CONT_motus_fdrs.csv')

    fdrs = [s0_2_fdr, s0_3_fdr, s0_4_fdr]
    dfs = [s0_2, s0_3, s0_4]
    res = []
    for df,name, fdr in zip(dfs, ['2', '3', '4'],fdrs ):
        temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
        temp.rename(columns={'Estimate': f'Estimate_{name}', 'Pr(>|t|)': f'Pr(>|t|)_{name}'}, inplace=True)
        fdr.columns = ['otu', f'q_{name}']
        temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
        temp['stage_vs_0'] = name
        res.append(temp)

    combo = res[0].merge(res[1], left_on='OTU', right_on='OTU', how='outer').merge(res[2], left_on='OTU', right_on='OTU', how='outer')
    combo = combo[['OTU', 'Estimate_2','Estimate_3', 'Estimate_4', 'Pr(>|t|)_2', 'Pr(>|t|)_3', 'Pr(>|t|)_4', 'q_2', 'q_3', 'q_4' ]].set_index('OTU')

    combo['if_2'] = combo['q_2'] < 0.1
    combo['if_3'] = combo['q_3'] < 0.1
    combo['if_4'] = combo['q_4'] < 0.1

    sub_combo = combo[(combo['if_2']==True) | (combo['if_3'] ==True) | (combo['if_4'] == True)]
    sub_combo['q_2'].values[sub_combo['q_2'] > 0.1] = np.nan
    sub_combo['q_3'].values[sub_combo['q_3'] > 0.1] = np.nan
    sub_combo['q_4'].values[sub_combo['q_4'] > 0.1] = np.nan

    sub_combo['log2'] = np.log10(sub_combo['q_2'])
    sub_combo['log3'] = np.log10(sub_combo['q_3'])
    sub_combo['log4'] = np.log10(sub_combo['q_4'])

    sub_combo['Stage 2 factor'] = sub_combo['log2'] * np.sign(sub_combo['Estimate_2'])
    sub_combo['Stage 3 factor'] = sub_combo['log3'] * np.sign(sub_combo['Estimate_3'])
    sub_combo['Stage 4 factor'] = sub_combo['log4'] * np.sign(sub_combo['Estimate_4'])

    motus_taxa = read_csv('Data/random/motus_taxa.csv')
    motus_taxa['otu name'] = motus_taxa['taxa'] + ' (' + motus_taxa['otu'] + ')'
    sub_combo = sub_combo.merge(motus_taxa, left_index=True, right_on='otu').set_index('otu name')
    motus_sub_combo = sub_combo.copy()
    motus_sub_combo.index = motus_sub_combo.index.str.replace('_', ' ').str.replace(' mOTU v25 ','_')
    motus_sub_combo.index = motus_sub_combo.index.str.replace('species incertae sedis ', 'sp. ')



    #PHAGES
    s0_2 = read_csv('Data/corncob/stage2_vs_CONT_phages.csv')
    s0_3 = read_csv('Data/corncob/stage3_vs_CONT_phages.csv')
    s0_4 = read_csv('Data/corncob/stage4_vs_CONT_phages.csv')

    s0_2_fdr = read_csv('Data/corncob/stage2_vs_CONT_phages_fdrs.csv')
    s0_3_fdr = read_csv('Data/corncob/stage3_vs_CONT_phages_fdrs.csv')
    s0_4_fdr = read_csv('Data/corncob/stage4_vs_CONT_phages_fdrs.csv')


    fdrs = [s0_2_fdr, s0_3_fdr, s0_4_fdr]
    dfs = [s0_2, s0_3, s0_4]
    res = []
    for df,name, fdr in zip(dfs, ['2', '3', '4'],fdrs ):
        temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
        temp.rename(columns={'Estimate': f'Estimate_{name}', 'Pr(>|t|)': f'Pr(>|t|)_{name}'}, inplace=True)
        fdr.columns = ['otu', f'q_{name}']
        temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
        temp['stage_vs_0'] = name
        res.append(temp)

    combo = res[0].merge(res[1], left_on='OTU', right_on='OTU', how='outer').merge(res[2], left_on='OTU', right_on='OTU', how='outer')
    combo = combo[['OTU', 'Estimate_2','Estimate_3', 'Estimate_4', 'Pr(>|t|)_2', 'Pr(>|t|)_3', 'Pr(>|t|)_4', 'q_2', 'q_3', 'q_4' ]].set_index('OTU')

    combo['if_2'] = combo['q_2'] < 0.1
    combo['if_3'] = combo['q_3'] < 0.1
    combo['if_4'] = combo['q_4'] < 0.1

    sub_combo = combo[(combo['if_2']==True) | (combo['if_3'] ==True) | (combo['if_4'] == True)]
    sub_combo['q_2'].values[sub_combo['q_2'] > 0.1] = np.nan
    sub_combo['q_3'].values[sub_combo['q_3'] > 0.1] = np.nan
    sub_combo['q_4'].values[sub_combo['q_4'] > 0.1] = np.nan

    sub_combo['log2'] = np.log10(sub_combo['q_2'])
    sub_combo['log3'] = np.log10(sub_combo['q_3'])
    sub_combo['log4'] = np.log10(sub_combo['q_4'])


    sub_combo['Stage 2 factor'] = sub_combo['log2'] * np.sign(sub_combo['Estimate_2'])
    sub_combo['Stage 3 factor'] = sub_combo['log3'] * np.sign(sub_combo['Estimate_3'])
    sub_combo['Stage 4 factor'] = sub_combo['log4'] * np.sign(sub_combo['Estimate_4'])

    phage_sub_combo = sub_combo.copy()

    ####PLASMIDS
    s0_2 = read_csv('Data/corncob/stage2_vs_CONT_plasmids.csv')
    s0_3 = read_csv('Data/corncob/stage3_vs_CONT_plasmids.csv')
    s0_4 = read_csv('Data/corncob/stage4_vs_CONT_plasmids.csv')

    s0_2_fdr = read_csv('Data/corncob/stage2_vs_CONT_plasmids_fdrs.csv')
    s0_3_fdr = read_csv('Data/corncob/stage3_vs_CONT_plasmids_fdrs.csv')
    s0_4_fdr = read_csv('Data/corncob/stage4_vs_CONT_plasmids_fdrs.csv')


    fdrs = [s0_2_fdr, s0_3_fdr, s0_4_fdr]
    dfs = [s0_2, s0_3, s0_4]
    res = []
    for df,name, fdr in zip(dfs, ['2', '3', '4'],fdrs ):
        temp = df.loc[df['Unnamed: 0'].str.startswith('mu.G')]
        temp.rename(columns={'Estimate': f'Estimate_{name}', 'Pr(>|t|)': f'Pr(>|t|)_{name}'}, inplace=True)
        fdr.columns = ['otu', f'q_{name}']
        temp = temp.merge(fdr, left_on='OTU', right_on='otu', how='left')
        temp['stage_vs_0'] = name
        res.append(temp)

    combo = res[0].merge(res[1], left_on='OTU', right_on='OTU', how='outer').merge(res[2], left_on='OTU', right_on='OTU', how='outer')
    combo = combo[['OTU', 'Estimate_2','Estimate_3', 'Estimate_4', 'Pr(>|t|)_2', 'Pr(>|t|)_3', 'Pr(>|t|)_4', 'q_2', 'q_3', 'q_4' ]].set_index('OTU')

    combo['if_2'] = combo['q_2'] < 0.1
    combo['if_3'] = combo['q_3'] < 0.1
    combo['if_4'] = combo['q_4'] < 0.1

    sub_combo = combo[(combo['if_2']==True) | (combo['if_3'] ==True) | (combo['if_4'] == True)]
    sub_combo['q_2'].values[sub_combo['q_2'] > 0.1] = np.nan
    sub_combo['q_3'].values[sub_combo['q_3'] > 0.1] = np.nan
    sub_combo['q_4'].values[sub_combo['q_4'] > 0.1] = np.nan

    sub_combo['log2'] = np.log10(sub_combo['q_2'])
    sub_combo['log3'] = np.log10(sub_combo['q_3'])
    sub_combo['log4'] = np.log10(sub_combo['q_4'])


    sub_combo['Stage 2 factor'] = sub_combo['log2'] * np.sign(sub_combo['Estimate_2'])
    sub_combo['Stage 3 factor'] = sub_combo['log3'] * np.sign(sub_combo['Estimate_3'])
    sub_combo['Stage 4 factor'] = sub_combo['log4'] * np.sign(sub_combo['Estimate_4'])


    plasmid_sub_combo = sub_combo.copy()

    sns.set(font='Seravek')
    plt.figure(figsize=(4,8))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    fig, axes = plt.subplots(3, 1,
            sharex=True,
            figsize=(1,8),
            gridspec_kw={'height_ratios': [len(motus_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0)),
                                           len(phage_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0)),
                                           len(plasmid_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0))]})
    sns.heatmap(
            motus_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0).sort_values(by=['Stage 2 factor','Stage 3 factor','Stage 4 factor']),
                       center=0, cmap=cmap, vmax=4,
            ax=axes[0],
            linewidths=0.01,
            yticklabels=False,
            cbar=False)
    sns.heatmap(
            phage_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0).sort_values(by=['Stage 2 factor','Stage 3 factor','Stage 4 factor']),
                       center=0, cmap=cmap, vmax=4,
            ax=axes[1],
            linewidths=0.01,
            yticklabels=False,
            cbar=False)

    sns.heatmap(
            plasmid_sub_combo[['Stage 2 factor','Stage 3 factor','Stage 4 factor']].fillna(0).sort_values(by=['Stage 2 factor','Stage 3 factor','Stage 4 factor']),
                       center=0, cmap=cmap, vmax=4,
            ax=axes[2],
            linewidths=0.01,
            yticklabels=False,
            cbar=False)


    axes[0].set_ylabel("")
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")
    axes[2].set_xlabel("")

    plt.savefig('/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/Fig4G.pdf', dpi=300)



if __name__ == "__main__":

    main()
