from numpy import unique
import flowkit as fk
import phenograph
import fcsparser
import xfcs
import glob
import os
import pandas as pd
from numpy import genfromtxt
import optparse
import sys
#
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-f', action="store", dest="input_filename", help='')
parser.add_option('-k', action="store", dest="kcoef", help='')
parser.add_option('-t', action="store", dest="thread", help='')
options, args = parser.parse_args()

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

def runphenograph(inputmatrix, outputfolder,prefix,kcoev,thread):
    """

    :param inputmatrix:
    :param outputfolder:
    :param prefix:
    :param kcoev:
    :param thread:
    :return:
    """
    # read input matrix
    inputmatrixpheno = pd.read_csv(inputmatrix,header=0)
    # run phenograph clustering
    communities, graph, Q = phenograph.cluster(inputmatrixpheno[inputmatrixpheno.columns.difference(['filename_cell'])], k=int(kcoev),
                                               directed=False,
                                               n_jobs=thread)
    # create dataframe with cluster information
    dfPheno = pd.DataFrame(communities)
    # rename colum 0
    dfPheno = dfPheno.rename(columns={0: 'cluster'})
    # merge input dataframe with phenograph clusters
    dfMerge = inputmatrixpheno.merge(dfPheno, how='outer', left_index=True, right_index=True)
    # export output
    dfMerge.to_csv("".join([outputfolder,prefix,"_clustered_data.csv"]),sep="\t",header=True)
    print('Found {} clusters'.format(len(unique(communities))), flush=True)
    # create a df groupby with cluster information
    cluster_frame = dfMerge.groupby("cluster", as_index=False).median()
    # skip -1 which means under min cluster size
    cluster_frame = cluster_frame[cluster_frame["cluster"] != -1] #
    cluster_frame["cluster"] = dfMerge.groupby("cluster")["cluster"].count()
    cluster_frame.to_csv("".join([outputfolder,prefix,"_cluster_info.csv"]),sep="\t",header=True)
    print('Clustering successful', flush=True)


if __name__ == '__main__':
    print("Clustering start")
    inputpath = pathchecker(options.input_folder)
    outpath = pathchecker(options.output_folder)
    runphenograph(inputmatrix="".join([inputpath,options.input_filename]),
                  outputfolder=outpath,
                  prefix=options.input_filename.split(".")[0],
                  kcoev=75,
                  thread=20)