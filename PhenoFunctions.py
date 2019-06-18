import glob
import os
import sys
import pandas as pd
import phenograph as pg
from openTSNE import TSNE
import flowkit as fk


def loadmarkers(markerfile):
    """
    Read marker filename with path
    :param markerfile:
    :return: array with filepath and filename
    """
    if (os.path.exists(markerfile) is True):
        markerfilepath = os.path.split(markerfile)[0]
        markerfilename = os.path.split(markerfile)[1]
        return (markerfilepath, markerfilename)
    else:
        print("File does not exist. Exit")
        sys.exit(1)


def loadconcatenatefile(concfile):
    """
    Read concatenate file with path
    :param concfile: output of inputaparsin.py (all markers)
    :return: array with filepath and filename
    """
    if (os.path.exists(concfile) is True):
        concfilepath = os.path.split(concfile)[0]
        concfilename = os.path.split(concfile)[1]
        return (concfilepath, concfilename)
    else:
        print("File does not exist. Exit")
        sys.exit(1)

def checkmarkers(markerfile,data):
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


def fcs2csv(inputfcsfolder,outputcsvfolder):
    """

    :param inputfcsfolder:
    :param outputcsvfolder:
    :return:
    """
    # set the extension value
    extension = 'fcs'
    # save path and input filename
    allfilenames = [i for i in glob.glob("".join([inputfcsfolder, '*.{}'.format(extension)]))]
    # save prefix name
    prefixname = []
    # convert with flowkit
    for i in range(len(allfilenames)):
        sample = fk.Sample(allfilenames[i], subsample_count=None)
        sample.export_csv(source='raw', subsample=False,
                          filename="".join([allfilenames[i].split("/")[-1].split(".")[0], ".csv"]),
                          directory=outputcsvfolder)


def csv2fcs(inputdf):
    """

    :param inputdf:
    :return:
    """
    sample = fk.Sample(inputdf.values, channel_labels=inputdf.columns)
    sample.export_fcs(source='raw',
                      filename="test_All.fcs",
                      directory="/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Jasper/CD4_mincsize1_lowersizeevent25/30/")


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
        df = df.append(file_df)
    df['event'] = df.index + 1
    df = df.reset_index(drop=True)
    # subsetting
    dfFileNameEvent = df[['file_name', 'event']]
    dfConcatenate = df.drop(columns=['file_name', 'event'])
    dfConcwitName = pd.merge(dfFileNameEvent, dfConcatenate,
                             left_index=True,
                             right_index=True, how='left')
    return dfFileNameEvent, dfConcatenate,dfConcwitName


def runphenograph(marker,concdf,kcoev,thread,outfold,name):
    """

    :param marker:
    :param condf:
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
    communities, graph, Q = pg.cluster(data.values, k=int(kcoev),
                                       directed=False,
                                       min_cluster_size=1,
                                       n_jobs=thread)
    dfPheno = pd.DataFrame(communities)
    dfMerge = pd.merge(data, dfPheno, left_index=True, right_index=True,    how = 'left')
    dfPheno.to_csv("/".join([outfold, "".join([name, "_clusters.txt"])]),
                   sep="\t",
                   header=True, index=False)
    dfMerge = dfMerge.rename(columns={0: 'Phenograph'})
    #dfMerge.columns = ['\"' + str(col) + '\"' for col in dfMerge.columns]
    dfMerge.to_csv("/".join([outfold, "".join([name, "_pheno.txt"])]),
                   sep="\t",
                   header=True, index=False)
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
                early_exaggeration_iter=250, early_exaggeration=12, n_iter=750,
                neighbors="exact", negative_gradient_method="bh")
    embedding = tsne.fit(x)
    dftsne = pd.DataFrame({'Tsne_1': embedding[:, 0], 'Tsne_2': embedding[:, 1]})
    dftsne.to_csv("/".join([outfold, "".join([name, "_tsne2.txt"])]),
                  sep="\t",
                  header=True, index=False)
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
    dfAll.to_csv("/".join([outfold, "".join([name, "_all.txt"])]),
                 sep=",",
                 header=True, index=False)
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
    df = df.drop([0, 1])
    head = df[df.columns[-1]].unique().tolist()
    for i in head:
        df.loc[df[df.columns[-1]] == i][df.loc[df[df.columns[-1]] == i].columns.difference([0,1])].to_csv("/".join([outfold, "FCScluster", "".join([name,"_",str(i),".csv"])]),
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
    for i in head:
        df.loc[df[df.columns[0]] == i][df.loc[df[df.columns[0]] == i].columns.difference([-1, -2])].to_csv("/".join([outfold, "FCSsample","".join([str(i),"_", name, ".csv"])]),
                                                                                                           header=True,
                                                                                                           index=False)

