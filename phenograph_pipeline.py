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
parser.add_option('-c', action="store", dest="channel_excluded", help='')
parser.add_option('-f', type='choice', choices=['csv', 'fcs', 'FCS', 'CSV'], help='')
options, args = parser.parse_args()


#
CSVraw_folder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/*.csv"
FCSraw_folder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs"
fnames = glob.glob(FCSraw_folder)
#print(fnames)


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


def read_markers(markerfile):
    """
    read file with marker to exclude and save in array
    :param markerfile:
    :return:
    """
    with open(markerfile) as my_file:
        marker_array = my_file.readlines()
    return marker_array


def csvmarkerremove(dataframeconcatenate,markerarray):

    #df = dataframeconcatenate
    for i in range(len(markerarray)):
        dataframeconcatenate.drop(markerarray[i].strip('\n'), axis=1, inplace=True)
    return dataframeconcatenate


def csvconcatenate(csvfolder):
    """

    :param csvfolder:
    :return:
    """
    extension = 'csv'
    os.chdir(csvfolder)
    allfilenames = [i for i in glob.glob('*.{}'.format(extension))]
    names = [os.path.basename(x) for x in glob.glob(csvfolder + '\*.csv')]
    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f, header=0) for f in allfilenames])
    # export to csv
    #combined_csv.to_csv("combined_csv.csv", index=False, header=True)

    all_files = glob.glob(os.path.join(csvfolder, "*.csv"))
    names = [os.path.basename(x) for x in all_files]
    df = pd.DataFrame()

    for file_, name in zip(all_files, names):
        file_df = pd.read_csv(file_, header=0)
        file_df['file_name'] = name.split('.')[0]
        df = df.append(file_df)

    df['index'] = df.index + 1
    df = df.set_index('index')
    df['index'] = df['file_name'] + '_' + df.index.astype(str)
    #df = df.set_index('index')
    df = df.drop(['file_name'], axis=1)

    return df


def csv2fcsconverter(csvfile, outputfolder):
    """

    :param csvfile:
    :param outputfolder:
    :return:
    """
    my_data = genfromtxt(csvfile, delimiter=',', skip_header=True)
    df = pd.read_csv(csvfile, header=0)
    sample = fk.Sample(fcs_path_or_data=my_data, subsample_count=None, channel_labels=df.columns)
    #print(sample)
    #print(sample.pnn_labels)
    fk.Sample.export_fcs(sample,
                         source='raw',
                         subsample=False,
                         filename=".".join([csvfile.split("/")[-1].split(".")[0], "fcs"]),
                         directory="".join([outputfolder, "/"]))
    return "/".join([outputfolder, ".".join([csvfile.split("/")[-1].split(".")[0], "fcs"])])


def runphenograph(inputmatrixpheno):
    communities, graph, Q = phenograph.cluster(inputmatrixpheno[inputmatrixpheno.columns.difference(['index'])], k=75,
                                               directed=False,
                                               n_jobs=32)
    print(pd.DataFrame(communities))
    dfPheno = pd.DataFrame(communities)
    dfPheno['index'] = dfPheno.index + 1
    dfPheno = dfPheno.rename(columns={'0': 'cluster'})
    dfPheno.to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/dfPheno.csv",sep="\t",header=True)
    inputmatrixpheno.to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/inputmatrixpheno.csv",sep="\t",header=True)
    # df_merged = inputmatrixpheno.merge(dfPheno, how='outer', left_index=True, right_index=True)
    # pd.DataFrame(communities).to_csv(
    #     "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/communities.csv")
    # df_merged.to_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/afterpheno.csv",
    #                  sep="\t",header=True,index=True)


#csv2fcsconverter(testcsv,FCSraw_folder)


if __name__ == '__main__':
    print("Script start")
    raw_folder = pathchecker(options.input_folder)
    output_folder = pathchecker(options.output_folder)
    createdir(FCSraw_folder)
    if options.channel_excluded is None:
        removestep = False
    else:
        removestep = True
        marker = read_markers(options.channel_excluded)
    print(marker)
    if options.f.upper() == "CSV":
        print("csv")
        if removestep is False:
            print("remove step false")
        else:
            pg_input = csvmarkerremove(csvconcatenate(raw_folder), marker)
            print(pg_input)
            runphenograph(pg_input)
    elif options.f.upper() == "FCS":
        print("fcs")
    else:
        print("Please specify input format: csv or fcs")
        sys.exit(1)


