import glob
import os
import sys
import pandas as pd
#import phenograph as pg
from openTSNE import TSNE
import flowkit as fk
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.cm as cm
import numpy as np
from sklearn.metrics import silhouette_samples, silhouette_score
import random
import subprocess
from numpy import unique
import shutil
random.seed(123456)
np.random.seed(123456)

sns.set_style('darkgrid')
sns.set_palette('muted')
sns.set_context("notebook", font_scale=1.5,
                rc={"lines.linewidth": 2.5})
RS = 123


def loadmarkers(markerfile):
    """
    Read marker filename with path
    :param markerfile:
    :return: array with filepath and filename
    """
    if os.path.exists(markerfile) is True:
        markerfilepath = os.path.split(markerfile)[0]
        markerfilename = os.path.split(markerfile)[1]
        return markerfilepath, markerfilename
    else:
        print("File does not exist. Exit")
        sys.exit(1)


def loadconcatenatefile(concfile):
    """
    Read concatenate file with path
    :param concfile: output of inputaparsin.py (all markers)
    :return: array with filepath and filename
    """
    if os.path.exists(concfile) is True:
        concfilepath = os.path.split(concfile)[0]
        concfilename = os.path.split(concfile)[1]
        return concfilepath, concfilename
    else:
        print("File does not exist. Exit")
        sys.exit(1)


def checkmarkers(markerfile, data):
    """
    Check if marker in file is also a column of conc file
    :param markerfile:
    :param data:
    :return:
    """
    # read marker file
    marker_array = [line.rstrip() for line in open(markerfile)]
    # read concatenate file
    for i in range(len(marker_array)):
        if {marker_array[i]}.issubset(data.columns):
            continue
        else:
            print("Marker {} not found in Matrix.".format(marker_array[i]))
            sys.exit(1)

    return marker_array


def fcs2csv(inputfcsfolder, outputcsvfolder):
    """

    :param inputfcsfolder:
    :param outputcsvfolder:
    :return:
    """
    # set the extension value
    extension = 'fcs'
    # save path and input filename
    allfilenames = [i for i in glob.glob("".join([inputfcsfolder, '*.{}'.format(extension)]))]
    # convert with flowkit
    for i in range(len(allfilenames)):
        sample = fk.Sample(allfilenames[i], subsample_count=None)
        sample.export_csv(source='raw', subsample=False,
                          filename="".join([allfilenames[i].split("/")[-1].split(".")[0], ".csv"]),
                          directory=outputcsvfolder)


def csv2fcs(inputdf, outfold, name):
    """

    :param inputdf:
    :param outfold:
    :param name:
    :return:
    """
    sample = fk.Sample(inputdf.values, channel_labels=inputdf.columns)
    sample.export_fcs(source='raw',
                      filename="".join([name, "_concatenated.fcs"]),
                      directory=outfold)


def filenameeventinfor(csvfolder):
    """

    :param csvfolder:
    :return:
    """
    # change folder
    os.chdir(csvfolder)
    # create array with all csvfile name and path
    all_files = glob.glob(os.path.join(csvfolder, "*.csv"))
    # get files name
    names = [os.path.basename(x) for x in all_files]
    # open a empty dataframe
    df = pd.DataFrame()
    # append csv files and add a new columns with file name of provenance
    for file_, name in zip(all_files, names):
        file_df = pd.read_csv(file_, header=0)
        file_df['file_name'] = name.split('.')[0]
        df = df.append(file_df,sort=False)
    df['event'] = df.index + 1
    df = df.reset_index(drop=True)
    # subsetting
    dfFileNameEvent = df[['file_name', 'event']]
    dfConcatenate = df.drop(columns=['file_name', 'event']) #TODO check col name
    dfConcwitName = pd.merge(dfFileNameEvent, dfConcatenate,
                             left_index=True,
                             right_index=True, how='left')
    return dfFileNameEvent, dfConcatenate, dfConcwitName


def runphenograph(marker, concdf, kcoev, thread, outfold, name):
    """

    :param marker:
    :param concdf:
    :param kcoev:
    :param thread:
    :param outfold:
    :param name:
    :return:
    """
    print("Start Phenograph")
    data = concdf.copy()
    for i in range(len(marker)):
        data.drop(marker[i], axis=1, inplace=True)
    data.to_csv("/".join([outfold, "".join([name, "_clusters.txt"])]), sep="\t",
                header=True, index=False)
    Rscript = shutil.which('Rscript')
    PhenoR = "/".join([os.path.dirname(os.path.abspath(__file__)), "Phenograph.R"])
    subprocess.call(['{0} {1} -i {2} -k {3} -o {4} -n {5}'.format(Rscript, PhenoR,
                                                                  "/".join([outfold, "".join([name, "_clusters.txt"])]),
                                                                  kcoev, outfold, name)], shell=True)
    # communities, graph, Q = pg.cluster(data.values, k=int(kcoev),
    #                                    directed=True,
    #                                    min_cluster_size=1,
    #                                    n_jobs=thread)
    # dfPheno = pd.DataFrame(communities)
    # dfMerge = pd.merge(data, dfPheno, left_index=True, right_index=True, how='left')
    # dfPheno.to_csv("/".join([outfold, "".join([name, "_clusters.txt"])]),
    #                sep="\t",
    #                header=True, index=False)
    # dfMerge = dfMerge.rename(columns={0: 'Phenograph'})
    # dfMerge.to_csv("/".join([outfold, "".join([name, "_pheno.txt"])]),
    #                sep="\t",
    #                header=True, index=False)
    dfPheno = pd.read_csv("/".join([outfold,"".join([name,".tsv"])]),sep="\t",header=0)
    # print('Found {} clusters'.format(len(unique(dfPheno["Phenograph"]))), flush=True)
    # create a df groupby with cluster information
    cluster_frame = dfPheno.groupby("Phenograph", as_index=False).median()
    cluster_frame.index += 1
    # skip -1 which means under min cluster size
    cluster_frame = cluster_frame[cluster_frame["Phenograph"] != -1]  #
    cluster_frame["Phenograph"] = dfPheno.groupby("Phenograph")["Phenograph"].count()
    cluster_frame.to_csv("".join([outfold, name, "_cluster_info.csv"]), sep="\t", header=True)
    print('', flush=True)
    print('Clustering successful', flush=True)
    return dfPheno


def runtsne(marker,concdf,thread,outfold,name):
    """

    :param marker:
    :param concdf:
    :param thread:
    :param outfold:
    :param name:
    :return:
    """
    print("Start tSNE")
    data = concdf.copy()
    for i in range(len(marker)):
        data.drop(marker[i], axis=1, inplace=True)
    x = data.values
    tsne = TSNE(n_components=2, perplexity=30, learning_rate=200,
                n_jobs=int(thread), initialization="pca", metric="euclidean",
                early_exaggeration_iter=250, early_exaggeration=12, n_iter=1000,
                neighbors="exact", negative_gradient_method="bh")
    embedding = tsne.fit(x)
    dftsne = pd.DataFrame({'Tsne_1': embedding[:, 0], 'Tsne_2': embedding[:, 1]})
    dftsne = dftsne.round(2)
    dfMerge = pd.merge(data, dftsne, left_index=True, right_index=True,
                       how='left')
    dfMerge.to_csv("/".join([outfold, "".join([name, "_tsne.txt"])]),
                   sep="\t",
                   header=True, index=False)
    return dftsne


def combine(data, tsne, pheno, outfold, name):
    """

    :param data:
    :param tsne:
    :param pheno:
    :param outfold:
    :param name:
    :return:
    """
    dftsnepheno = pd.merge(tsne, pheno, left_index=True, right_index=True,
                           how='left')
    dfAll = pd.merge(data, dftsnepheno, left_index=True, right_index=True,
                     how='left')
    dfAll = dfAll.rename(columns={0: 'Phenograph'})
    dfAll.columns = ['\"' + str(col) + '\"' for col in dfAll.columns]
    dfout = dfAll.copy()
    removec = []
    removec += [col for col in dfout.columns if 'event' in col]
    removec += [col for col in dfout.columns if 'file_name' in col]
    dfout = dfout.drop(columns=[removec[0], removec[1]])
    dfout.to_csv("/".join([outfold, "".join([name, "_concatenated.txt"])]),
                 sep=",",
                 header=True, index=False)
    csv2fcs(inputdf=dfout, outfold=outfold, name=name)
    return dfAll


def createdir(dirpath):
    """
    Make dir function and check if directory is already exists
    :param dirpath: string with path and directory name
    :return:
    """
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
    else:
        print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))


def groupbycluster(alldf, outfold, name):
    """

    :param alldf:
    :param outfold:
    :param name:
    :return:
    """
    createdir("/".join([outfold, "FCScluster"]))
    df = alldf.copy()
    head = df[df.columns[-1]].unique().tolist()
    removec = []
    removec += [col for col in df.columns if 'event' in col]
    removec += [col for col in df.columns if 'file_name' in col]
    for i in head:
        df.loc[df[df.columns[-1]] == i].drop(columns=[removec[0], removec[1]]).to_csv(
            "/".join([outfold, "FCScluster", "".join([name, "_", str(i), ".csv"])]),
            header=True, index=False)


def groupbysample(alldf, outfold, name):
    """

    :param alldf:
    :param outfold:
    :param name:
    :return:
    """
    createdir("/".join([outfold, "FCSsample"]))
    df = alldf.copy()
    head = df[df.columns[0]].unique().tolist()
    removec = []
    removec += [col for col in df.columns if 'event' in col]
    removec += [col for col in df.columns if 'file_name' in col]
    for i in head:
        df.loc[df[df.columns[0]] == i].drop(columns=[removec[0], removec[1]]).to_csv(
            "/".join([outfold, "FCSsample", "".join([str(i), "_", name, ".csv"])]),
            header=True, index=False)


def validationplot(marker,alldf,outfold,name):
    """

    :param marker:
    :param alldf:
    :param outfold:
    :param name:
    :return:
    """
    data = alldf.copy()
    tsne = alldf.copy()
    head = data[data.columns[0]].unique().tolist()
    removec = []
    removec += [col for col in data.columns if 'event' in col]
    removec += [col for col in data.columns if 'file_name' in col]
    removec += [col for col in data.columns if 'Phenograph' in col]
    removec += [col for col in data.columns if 'Tsne_1' in col]
    removec += [col for col in data.columns if 'Tsne_2' in col]
    data.drop(columns=[removec[0],removec[1]], axis=1, inplace=True)
    data = data.loc[data[removec[2]] >= 0]
    for i in range(len(marker)):
        data.drop("\""+marker[i]+"\"", axis=1, inplace=True)
    y = data[removec[2]].values.flatten()
    X = data.drop(columns=[removec[2],removec[3],removec[4]]).values
    # n_clusters = [int(data[removec[2]].max() + 1)][0] old command
    n_clusters = [int(data[removec[2]].max())][0]
    # fig, (ax1, ax2) = plt.subplots(1, 1)
    # fig.set_size_inches(18, 7)
    f = plt.figure(figsize=(18, 7))
    plt.style.use('seaborn-white')
    plt.xlim(-1, 1)
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    plt.ylim(0, len(X) + (n_clusters + 1) * 10)


    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    # clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = y
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(1,n_clusters+1):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        plt.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    #plt.title("The silhouette plot for the various clusters.")
    plt.xlabel("The silhouette coefficient values ")
    plt.ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    plt.axvline(x=silhouette_avg, color="red", linestyle="--")

    plt.yticks([])  # Clear the yaxis labels / ticks
    plt.xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    plt.title(("Silhouette analysis for Phenograph clustering on {0} data with n_clusters = {1}.\n The average silhouette_score is {2}".format(name,n_clusters,silhouette_avg)),
                 fontsize=14, fontweight='bold')
    plt.savefig("/".join([outfold, "".join([name, "_validation"])]))


def tsneplot(x, colors,outfold, name):
    """

    :param x:
    :param colors:
    :param outfold:
    :param name:
    :return:
    """
    # choose a color palette with seaborn.
    num_classes = len(np.unique(colors))
    palette = np.array(sns.color_palette("hls", num_classes+1))
    #print(colors.max())
    #print((colors.min()))
    #print(palette)
    # create a scatter plot.
    f = plt.figure(figsize=(8, 8))
    ax = plt.subplot(aspect='equal')
    scatter = ax.scatter(x[:,0], x[:,1], lw=0, s=40,alpha=0.3, c=palette[colors.values.astype(np.int)])
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)
    ax.axis('off')
    ax.axis('tight')
    # add the labels for each digit corresponding to the label
    txts = []
    for i in range(1,num_classes+1):

        # Position of each label at median of data points.

        xtext, ytext = np.median(x[colors == i, :], axis=0)
        txt = ax.text(xtext, ytext, str(i), fontsize=24)
        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        txts.append(txt)

    # produce a legend with the unique colors from the scatter
    plt.savefig("/".join([outfold, name]))
    return f, ax, scatter, txts
