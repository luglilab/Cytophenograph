import numpy as np
from openTSNE import TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import optparse
plt.switch_backend('agg')
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set(style="white")
#
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-c', action="store", dest="channel_excluded", help='')
parser.add_option('-n', action="store", dest="analysis_name", help='')
parser.add_option('-t', action="store", dest="thread", help='')
options, args = parser.parse_args()
#
def pathchecker(inputpath):
    """
    Fix input or output string output
    :param inputpath: path to check
    :return: string
    """
    if inputpath.endswith('/') is True:
        return inputpath
    else:
        return "".join([inputpath, "/"])


def tsneplot(inpfile,thread, outpath, name):
    """

    :param inpfile:
    :param thread:
    :param outpath:
    :param name:
    :return:
    """
    dfOutPheno = pd.read_csv(inpfile, header=0, sep="\t", index_col="filename_cell")
    dfOutPheno = dfOutPheno.drop(['Unnamed: 0'], axis=1)
    dfTarget = dfOutPheno[['cluster']]
    x = dfOutPheno.values
    y = dfTarget.values
    tsne = TSNE(n_components=2, perplexity=30, learning_rate=200,
                n_jobs=int(thread), initialization="pca", metric="euclidean",
                early_exaggeration_iter=250, early_exaggeration=12, n_iter=750,
                neighbors="exact", negative_gradient_method="bh")
    embedding = tsne.fit(x)
    print(embedding)
    # # np.savetxt("foo.csv", embedding, delimiter=",")
    # dftsne = pd.DataFrame({'Tsne_1': embedding[:, 0], 'Tsne_2': embedding[:, 1]})
    # dfOutPheno = pd.read_csv(inpfile, header=0, sep="\t", index_col="Unnamed: 0")
    # # print(dfOutPheno.head())
    # # print(dataset.head())
    # dfMerge = pd.merge(dfOutPheno, dftsne, left_index=True, right_index=True, how = 'left')
    # # print(dfMerge.head())
    # # df = pd.read_csv("foo.csv", sep=",", header=None)
    # #
    # # print(df.head())
    # # tips = sns.load_dataset("tips")
    # # print(tips.head())
    # # dfMerge.to_csv("tsneinput.txt", sep="\t", header=True, index=False)
    # # dfMerge = pd.read_csv("tsneinput.txt",sep="\t",header=0)
    # # pal = sns.palplot(sns.color_palette("muted"))
    # testplot = sns.relplot(x="Tsne_1", y="Tsne_2", hue="cluster", palette='Spectral',
    #                        data=dfMerge, legend="full", sizes=10)
    # testplot.savefig("".join([outpath, name]))


if __name__ == '__main__':
    #inputpath = pathchecker(options.input_folder)
    outpath = pathchecker(options.output_folder)
    tsneplot(inpfile=options.input_folder,thread=options.thread,outpath=outpath, name=options.analysis_name)







