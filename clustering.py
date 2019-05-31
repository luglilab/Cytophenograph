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


def runphenograph(inputmatrix, outputfolder):
    inputmatrixpheno = pd.read_csv(inputmatrix,header=0,index_col="filename_cell")
    print(inputmatrixpheno.head())
    communities, graph, Q = phenograph.cluster(inputmatrixpheno[inputmatrixpheno.columns.difference(['index','filename_cell'])], k=75,
                                               directed=False,
                                               n_jobs=32)
    dfPheno = pd.DataFrame(communities)
    dfPheno['index'] = dfPheno.index + 1
    dfPheno = dfPheno.rename(columns={'0': 'cluster'})
    #dfPheno.to_csv("".join([outputfolder,"dfPheno.csv"]),sep="\t",header=True)
    #inputmatrixpheno.to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/inputmatrixpheno.csv",sep="\t",header=True)
    #print(dfPheno.head())
    #print(inputmatrixpheno.head())
    print('Found {} clusters'.format(len(unique(communities))), flush=True)
    new = inputmatrixpheno.index.str.split("_", n=0, expand=True)
    inputmatrixpheno['event'] = new.iloc[:, -1]
    df_merged = pd.merge(inputmatrixpheno,dfPheno,left_on="event",right_on="index")
    pd.DataFrame(communities).to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/communities.csv")
    df_merged.to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/afterpheno.csv",
                     sep="\t",header=True,index=True)


if __name__ == '__main__':
    runphenograph("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs/test_nomarkerscaled.txt",
                  outputfolder="/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/")