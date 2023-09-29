###data from: "/Users/ab4966/Downloads/signif_vsvs_annotations_sumary.xlsx"
from pandas import DataFrame, melt
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chisquare
from os.path import join

outdir = '/Users/ab4966/Dropbox/2022_CNIO/Figures/Figure_parts/'


def main():
    plt.close('all')
    fig = plt.figure(figsize = (3,2))

     ####aggregated bytes SV
    res = DataFrame([[14,22],[7,0], [8,7]], columns=['assembled flanking','reference genomes'], index=['agreement', 'disagreement', 'un-annoted'])


    g=sns.barplot(x='index', y='value', data=melt(res.reset_index(), id_vars='index',value_vars=['assembled flanking','reference genomes']),
                    hue='variable', palette=['mediumseagreen','mediumseagreen'], width=0.7, edgecolor='black')

    plt.xticks(color='w')
    plt.yticks(color='w')
    plt.xlabel('')
    plt.ylabel('')
    plt.legend([],[], frameon=False)
    plt.savefig(join(outdir,'fig1E_bar.pdf'),bbox_inches='tight', dpi=600)
    print(chisquare(res, axis=None))


if __name__ == "__main__":

    main()

