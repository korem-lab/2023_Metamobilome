import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv
from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
from scipy import stats
import pandas as pd

def add_p_val(lft,rgt,y,h,p):
    plt.plot([lft, lft, rgt, rgt], [y, y+h, y+h, y], lw=1, c='k')
    plt.text((lft + rgt) * .5, y+h, ('n.s.' if p > 0.15 else 'p < %.2g' if p > 0.001 else 'p < %.1g') % max(p+1e-20, 1e-20), ha='center', va='bottom', color='k',fontsize=6, fontname='Seravek')

def read_meta_cnio(meta_file):
    meta = read_csv(meta_file, sep='\t')
    meta['source'] = meta.sample_alias.str.split('-').str[2]
    meta = meta.set_index('sample_alias')
    meta = meta[~meta.index.duplicated(keep='first')]
    meta = meta.loc[meta.subject_disease_status != 'NonElegible']
    meta = meta.dropna(subset=['diabin'])
    meta.diabin = meta.diabin.astype('int')
    return meta


meta = read_meta_cnio('/Users/ab4966/OneDrive - cumc.columbia.edu/Eukarotes_calling_from_metagenomes/CNIO_sub/ids_metaG_selected_Variables.tsv')


def compare_aucs(au1, au2):
    """
    au1- list of aucs
    au2 - other list of aucs
    """
    au1 = np.array(au1)
    au2 = np.array(au2)

    m1, m2 = au1.mean(), au2.mean()
    r = stats.pearsonr(au1, au2)[0]
    print(m1, m2, r)
    s1, s2 = np.std(au1), np.std(au2)

    z = (m1 - m2) / np.sqrt(np.power(s1, 2) + np.power(s2, 2) - r * s1 * s2)
    p = 2*stats.norm.sf(abs(z))
    print("p-value: " + str(p))
    print("------------------")
    return (p)



def plot_robust_rocs2_two_cohorts(md1, md2, dfs1, dfs2, plot_all_lines1, plot_all_lines2, labels1, labels2, colors1, colors2, figname):
    preds = {}
    cur_aucs = {}
    base_fpr = np.linspace(0, 1, 101)
    plt.close('all')
    plt.figure(figsize=(2.7, 2.7))

    for df, l, c, plot_l in zip(dfs1, labels1, colors1, plot_all_lines1):
        df = df.loc[md1.index]

        aucs = []
        tprs = []

        for col in df.columns:
            fpr, tpr, _ = roc_curve(md1['is_pc'], df[col])
            if plot_l:
                plt.plot(fpr, tpr, alpha=0.1, color=c)

            auc = roc_auc_score(md1['is_pc'], df[col])
            aucs.append(auc)

            tpr = np.interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)

        mean_auc = np.mean(aucs)
        df.rename(columns={col: l}, inplace=True)
        preds[l] = pd.concat([md1['is_pc'], df[l]], axis=1)
        cur_aucs[l] = aucs
        tprs = np.array(tprs)
        mean_tprs = tprs.mean(axis=0)
        plt.plot(base_fpr, mean_tprs, label=f"{l}, auROC={np.round(mean_auc, 2)}", color=c)

    for df, l, c, plot_l, in zip(dfs2, labels2, colors2, plot_all_lines2):
        df = df.loc[md2.index]

        aucs = []
        tprs = []

        for col in df.columns:
            fpr, tpr, _ = roc_curve(md2['is_pc'], df[col])
            if plot_l:
                plt.plot(fpr, tpr, alpha=0.1, color=c)

            auc = roc_auc_score(md2['is_pc'], df[col])
            aucs.append(auc)

            tpr = np.interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)

        mean_auc = np.mean(aucs)
        preds[l] = pd.concat([md2['is_pc'], df[col]], axis=1)
        cur_aucs[l] = aucs
        tprs = np.array(tprs)
        mean_tprs = tprs.mean(axis=0)
        plt.plot(base_fpr, mean_tprs, label=f"{l}, auROC={np.round(mean_auc, 2)}", color=c)


    plt.legend(fontsize=4)
    plt.plot([0, 1], [0, 1], '--', c='grey')
    plt.xticks(color='w')
    plt.yticks(color='w')

    plt.savefig(figname,bbox_inches='tight', dpi=600)

    return preds, cur_aucs
