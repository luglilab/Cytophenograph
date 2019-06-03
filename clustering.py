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





def runphenograph(inputmatrix, outputfolder,prefix):
    inputmatrixpheno = pd.read_csv(inputmatrix,header=0)
    communities, graph, Q = phenograph.cluster(inputmatrixpheno[inputmatrixpheno.columns.difference(['filename_cell'])], k=75,
                                               directed=False,
                                               n_jobs=20)
    dfPheno = pd.DataFrame(communities)
    dfPheno = dfPheno.rename(columns={0: 'cluster'})
    dfMerge = inputmatrixpheno.merge(dfPheno, how='outer', left_index=True, right_index=True)
    dfMerge.to_csv("".join([outputfolder,prefix,"_clustered_data.csv"]),sep="\t",header=True)
    print('Found {} clusters'.format(len(unique(communities))), flush=True)
    cluster_frame = dfMerge.groupby("cluster", as_index=False).median()
    cluster_frame = cluster_frame[cluster_frame["cluster"] != -1]
    cluster_frame["cluster"] = dfMerge.groupby("cluster")["cluster"].count()
    cluster_frame.to_csv("".join([outputfolder,prefix,"_cluster_info.csv"]),sep="\t",header=True)
    print('Clustering successful', flush=True)


if __name__ == '__main__':
    runphenograph("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs/test_nomarkerscaled.txt",
                  outputfolder="/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/",prefix="test1")