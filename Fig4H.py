import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from will_eval_functions import get_robust_preds
from os.path import join

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'

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



def plot_robust_rocs2(md, dfs, labels, colors, figname):

    base_fpr = np.linspace(0, 1, 101)
    plt.close('all')
    plt.figure(figsize=(2.7, 2.7))

    for df, l, c in zip(dfs, labels, colors):
        df = df.loc[md.index]

        aucs = []
        tprs = []

        for col in df.columns:
            fpr, tpr, _ = roc_curve(md['is_pc'], df[col])
            plt.plot(fpr, tpr, alpha=0.1, color=c)

            auc = roc_auc_score(md['is_pc'], df[col])
            aucs.append(auc)

            tpr = np.interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)

        mean_auc = np.mean(aucs)
        tprs = np.array(tprs)
        mean_tprs = tprs.mean(axis=0)
        plt.plot(base_fpr, mean_tprs, label=f"{l}, auROC={np.round(mean_auc, 2)}", color=c)


    plt.legend(fontsize=5)
    plt.plot([0, 1], [0, 1], '--', c='grey')
    plt.xticks(color='w')
    plt.yticks(color='w')

    plt.savefig(figname,bbox_inches='tight', dpi=600)


def main():

    custom_palette = sns.color_palette("Paired", 10)

    motus_raw = pd.read_pickle('Data/LightGBM_predictions/early_motus_best.pkl')
    plasmid_raw = pd.read_pickle('Data/LightGBM_predictions/early_plasmid_best_fixed_binary.pkl')
    phage_raw = pd.read_pickle('Data/LightGBM_predictions/early_phage_best.pkl')

    motus_preds = get_robust_preds(motus_raw[0][0], cv_n=10, classifier=True)
    plasmid_preds = get_robust_preds(plasmid_raw[0][0], cv_n=10, classifier=True)
    phage_preds = get_robust_preds(phage_raw[0][0], cv_n=10, classifier=True)

    s_md = pd.read_csv('Data/LightGBM_predictions/infiles/s_meta.csv')
    s_md = s_md.loc[~s_md.stage.isin([3.0, 4.0])]
    s_md = s_md.set_index('sample_id')
    dfs1 = [motus_preds, plasmid_preds, phage_preds]
    labels1 = ['Early-Bacteria-Spanish (CV)', 'Early-Plasmids-Spanish (CV)', 'Early-Phages-Spanish (CV)']
    colors1 = [custom_palette[1], custom_palette[3], custom_palette[7]]
    plotlines1 = [True, True, True]


    ###germans

    g_motus_tr_test = pd.read_pickle('Data/LightGBM_predictions/e_motus_tr_test.pkl')
    g_plasmid_tr_test = pd.read_pickle('Data/LightGBM_predictions/e_plas_tr_test.pkl')
    g_phage_tr_test = pd.read_pickle('Data/LightGBM_predictions/e_phag_tr_test.pkl')

    g_motus_preds_tt = get_robust_preds(g_motus_tr_test[0][0], cv_n=10, classifier=True)
    g_plasmid_preds_tt = get_robust_preds(g_plasmid_tr_test[0][0], cv_n=10, classifier=True)
    g_phage_preds_tt = get_robust_preds(g_phage_tr_test[0][0], cv_n=10, classifier=True)

    g_md = pd.read_csv('Data/LightGBM_predictions/infiles/german_meta.csv')
    g_md = g_md.loc[~g_md.stage.isin([3.0, 4.0])]
    g_md = g_md.set_index('sample_id')
    dfs2 = [g_motus_preds_tt, g_plasmid_preds_tt, g_phage_preds_tt]
    labels2 = ['Early Bacteria-German(Val.)', 'Early Plasmids-German(Val.)', 'Early Phages-German(Val.)']
    colors2 = [custom_palette[0], custom_palette[2], custom_palette[6]]
    plotlines2 = [True, True, True]

    cur_preds, cur_aucs = plot_robust_rocs2_two_cohorts(s_md, g_md, dfs1, dfs2, plotlines1, plotlines2, labels1, labels2, colors1, colors2, join(outdir,'/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/fig4H_auc.pdf'))


if __name__ == "__main__":
    main()
