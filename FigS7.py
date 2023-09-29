import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from will_eval_functions import eval_iterations, get_robust_preds


def plot_robust_rocs2_two_cohorts(md1, md2, dfs1, dfs2, plot_all_lines1, plot_all_lines2, labels1, labels2, colors1, colors2, figname):

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
        tprs = np.array(tprs)
        mean_tprs = tprs.mean(axis=0)
        plt.plot(base_fpr, mean_tprs, label=f"{l}, auROC={np.round(mean_auc, 2)}", color=c)


    plt.legend(fontsize=4)
    plt.plot([0, 1], [0, 1], '--', c='grey')
    plt.xticks(color='w')
    plt.yticks(color='w')

    plt.savefig(figname,bbox_inches='tight', dpi=600)

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

    s_md = pd.read_csv('Data/LightGBM_predictions/infiles/s_meta.csv')
    s_md = s_md.loc[~s_md.stage.isin([3.0, 4.0])]
    s_md = s_md.set_index('sample_id')

    g_md = pd.read_csv('Data/LightGBM_predictions/infiles/german_meta.csv')
    g_md = g_md.loc[~g_md.stage.isin([3.0, 4.0])]
    g_md = g_md.set_index('sample_id')


    motus_plas = pd.read_pickle('Data/LightGBM_predictions/early_motus_lra_plasmid_raw_best.pkl')
    motus_phag = pd.read_pickle('Data/LightGBM_predictions/early_motus_lra_phage_raw_best.pkl')
    trio = pd.read_pickle('Data/LightGBM_predictions/early_motus_lra_plasmid_raw_phage_raw_best.pkl')

    mts_plas_preds = get_robust_preds(motus_plas[0][0], cv_n=10, classifier=True)
    mts_phag_preds = get_robust_preds(motus_phag[0][0], cv_n=10, classifier=True)
    trio_preds = get_robust_preds(trio[0][0], cv_n=10, classifier=True)

    dfs1 = [mts_plas_preds, mts_phag_preds, trio_preds]
    labels1 = ['Early Bacteria and Plasmids(CV)', 'Early Bacteria and Phages(CV)', 'Early Bacteria,Plasmids and Phages(CV)']
    colors1 = [custom_palette[1], custom_palette[3], custom_palette[7]]
    plotlines1 = [True, True, True]

    #trained CNIO tested on GE

    motus_plas_tt = pd.read_pickle('Data/LightGBM_predictions/e_mts_plas_tr_test.pkl')
    motus_plas_preds_tt = get_robust_preds(motus_plas_tt[0][0], cv_n=10, classifier=True)

    motus_phag_tt = pd.read_pickle('Data/LightGBM_predictions/e_mts_phage_tr_test.pkl')
    motus_phage_preds_tt = get_robust_preds(motus_phag_tt[0][0], cv_n=10, classifier=True)

    trio_tt = pd.read_pickle('Data/LightGBM_predictions/e_mts_plas_phag_tr_test.pkl')
    trio_preds_tt = get_robust_preds(trio_tt[0][0], cv_n=10, classifier=True)

    dfs2 = [motus_plas_preds_tt, motus_phage_preds_tt, trio_preds_tt]
    labels2 = ['Early Bacteria and Plasmids-German(Val.)', 'Early Bacteria and Phages-German(Val.)', 'Early Bacteria,Plasmids and Phages-German(Val.)']
    colors2 = [custom_palette[0], custom_palette[2], custom_palette[6]]
    plotlines2 = [True, True, True]

    plot_robust_rocs2_two_cohorts(s_md, g_md, dfs1, dfs2, plotlines1, plotlines2, labels1, labels2, colors1, colors2, '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/figS7.pdf')


if __name__ == "__main__":
    main()
