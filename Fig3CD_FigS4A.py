import pandas as pd
import seaborn as sns
from will_eval_functions import get_robust_preds
from utils import plot_robust_rocs2_two_cohorts
from os.path import join

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'



def main():

    custom_palette = sns.color_palette("Paired", 10)

    motus_raw = pd.read_pickle('Data/LightGBM_predictions/s_motus_best.pkl')
    plasmid_raw = pd.read_pickle('Data/LightGBM_predictions/s_plasmid_best.pkl')
    phage_raw = pd.read_pickle('Data/LightGBM_predictions/s_phage_best.pkl')

    motus_preds = get_robust_preds(motus_raw[0][0], cv_n=10, classifier=True)
    plasmid_preds = get_robust_preds(plasmid_raw[0][0], cv_n=10, classifier=True)
    phage_preds = get_robust_preds(phage_raw[0][0], cv_n=10, classifier=True)

    clin_model = pd.read_pickle('Data/LightGBM_predictions/s_clin_data_hps.pkl')
    clin_preds = get_robust_preds(clin_model[0][0], cv_n=10, classifier=True)

    s_md = pd.read_csv('Data/LightGBM_predictions/infiles/s_meta.csv')
    s_md = s_md.set_index('sample_id')
    dfs1 = [motus_preds, plasmid_preds, phage_preds]
    labels1 = ['Bacteria-Spanish (CV)', 'Plasmids-Spanish (CV)', 'Phages-Spanish (CV)']
    colors1 = [custom_palette[1], custom_palette[3], custom_palette[7]]
    plotlines1 = [True,True,True]


    ###germans
    g_motus_tr_test = pd.read_pickle('Data/LightGBM_predictions/motus_tr_test.pkl')
    g_plasmid_tr_test = pd.read_pickle('Data/LightGBM_predictions/plasmid_tr_test.pkl')
    g_phage_tr_test = pd.read_pickle('Data/LightGBM_predictions/phage_tr_test.pkl')

    g_motus_preds_tt = get_robust_preds(g_motus_tr_test[0][0], cv_n=10, classifier=True)
    g_plasmid_preds_tt = get_robust_preds(g_plasmid_tr_test[0][0], cv_n=10, classifier=True)
    g_phage_preds_tt = get_robust_preds(g_phage_tr_test[0][0], cv_n=10, classifier=True)


    g_md = pd.read_csv('Data/LightGBM_predictions/infiles/german_meta.csv')
    g_md = g_md.set_index('sample_id')
    dfs2 = [g_motus_preds_tt, g_plasmid_preds_tt, g_phage_preds_tt]
    labels2 = ['Bacteria-German (Val.)', 'Plasmids-German (Val.)', 'Phages-German (Val.)']
    colors2 = [custom_palette[0], custom_palette[2], custom_palette[6]]
    plotlines2 = [True, True, True]

    cur_preds, cur_aucs = plot_robust_rocs2_two_cohorts(s_md, g_md, dfs1 +[clin_preds], dfs2, plotlines1 + [False], plotlines2,
                                                        labels1 + ['Clinical-Spanish (CV)'], labels2, colors1 + ['darkred'],
                                                        colors2, join(outdir,'fig3C_auc.pdf')
                                                        )

    ###Japanese indiv model
    j_motus_tr_test = pd.read_pickle('Data/LightGBM_predictions/motus_tr_jap_test.pkl')
    j_plasmid_tr_test = pd.read_pickle('Data/LightGBM_predictions/plasmid_tr_jap_test.pkl')
    j_phage_tr_test = pd.read_pickle('Data/LightGBM_predictions/phage_tr_jap_test.pkl')

    j_motus_preds_tt = get_robust_preds(j_motus_tr_test[0][0], cv_n=10, classifier=True)
    j_plasmid_preds_tt = get_robust_preds(j_plasmid_tr_test[0][0], cv_n=10, classifier=True)
    j_phage_preds_tt = get_robust_preds(j_phage_tr_test[0][0], cv_n=10, classifier=True)


    j_md = pd.read_csv('Data/LightGBM_predictions/infiles/japan_meta.csv')
    j_md = j_md.set_index('sample_id')
    dfs3 = [j_motus_preds_tt, j_plasmid_preds_tt, j_phage_preds_tt]
    labels3 = ['Bacteria-Japanese (Val.)', 'Plasmids-Japanese (Val.)', 'Phages-Japanese (Val.)']
    colors3 = [custom_palette[0], custom_palette[2], custom_palette[6]]
    plotlines3 = [True, True, True]

    cur_preds, cur_aucs = plot_robust_rocs2_two_cohorts(s_md, j_md, dfs1, dfs3, [False, False, False], plotlines3, labels1, labels3, colors1, colors3, join(outdir,'figS4A.pdf'))



################FIGURE 3D motus plus combos
    #spanish_combos
    motus_plas = pd.read_pickle('Data/LightGBM_predictions/motus_raw_plasmid_clr_best.pkl')
    motus_plas_preds = get_robust_preds(motus_plas[0][0], cv_n=10, classifier=True)

    trio = pd.read_pickle('Data/LightGBM_predictions/motus_raw_plasmid_clr_phage_lra_best.pkl')
    trio_pred = get_robust_preds(trio[0][0], cv_n=10, classifier=True)

    dfs1 = [motus_preds, motus_plas_preds, trio_pred]
    labels1 = ['Bacteria-Spanish(CV)', 'Bacteria and Plasmids-Spanish(CV)', 'Bacteria,Plasmid and Phages-Spanish(CV)']
    colors1 = [custom_palette[1], custom_palette[5], custom_palette[9]]
    plotlines1 = [False, True, True]

    #german_new

    motus_plas_tt = pd.read_pickle('Data/LightGBM_predictions/motus_plas_tr_test.pkl')
    motus_plas_preds_tt = get_robust_preds(motus_plas_tt[0][0], cv_n=10, classifier=True)

    trio_tt = pd.read_pickle('Data/LightGBM_predictions/motus_phag_plas_tr_test.pkl')
    trio_preds_tt = get_robust_preds(trio_tt[0][0], cv_n=10, classifier=True)



    dfs2 = [g_motus_preds_tt, motus_plas_preds_tt, trio_preds_tt]
    labels2 = ['Bacteria-German(Val.)', 'Bacteria and Plasmids-German(Val.)', 'Bacteria,Plasmid and Phages-German(Val.)']
    colors2 = [custom_palette[0], custom_palette[4], custom_palette[8]]

    plotlines2 = [False, True, True]

    cur_preds,cur_aucs = plot_robust_rocs2_two_cohorts(s_md, g_md, dfs1, dfs2, plotlines1, plotlines2, labels1, labels2, colors1, colors2, join(outdir, 'fig3D_auc.pdf'))


if __name__ == "__main__":
    main()
