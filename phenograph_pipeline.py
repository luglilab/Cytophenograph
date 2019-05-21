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
    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f,header=0) for f in allfilenames])
    # export to csv
    #combined_csv.to_csv("combined_csv.csv", index=False, header=True)
    return combined_csv


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
    print(options)
    print(raw_folder)
    print(output_folder)
    print(marker)
    if options.f.upper() == "CSV":
        print("csv")
        if removestep is False:
            print("remove step false")
        else:
            #csvconcatenate(raw_folder)
            pg_input = csvmarkerremove(csvconcatenate(raw_folder), marker)
            print(pg_input.head())
            communities, graph, Q = phenograph.cluster(pg_input, k=75, directed=False, n_jobs=32)
            print(communities)
            print(graph)
            print(Q)
    elif options.f.upper() == "FCS":
        print("fcs")
    else:
        print("Please specify input format: csv or fcs")
        sys.exit(1)


