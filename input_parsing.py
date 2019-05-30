import glob
import os
import pandas as pd
import optparse
import sys
import flowkit as fk
#
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-c', action="store", dest="channel_excluded", help='')
parser.add_option('-n', action="store", dest="analysis_name", help='')
parser.add_option('-f', type='choice', choices=['csv', 'fcs', 'FCS', 'CSV'], help='')
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
    Read file with marker to exclude and save in array
    :param markerfile:
    :return:
    """
    with open(markerfile) as my_file:
        marker_array = my_file.readlines()

    return marker_array

def csvconcatenate(csvfolder,outpfold,outputname):
    """
    Append all csv file in order to generate input matrix for phenograph
    :param csvfolder:
    :param outpfold:
    :return:
    """
    extension = 'csv'
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
    # create a new column called events with the cell number
    df['event'] = df.index + 1
    # reset the input
    df.reset_index(inplace=True)
    # index start from 1
    df.index += 1
    # create a column called filename_cell with file name and event
    df['filename_cell'] = df['file_name'] + '_' + df['event'].astype(str)
    # remove intermediate
    df = df.drop(['file_name', 'index'], axis=1)
    # export to csv
    df.to_csv("".join([outpfold,outputname,".csv"]), index=False, header=True)
    return "".join([outpfold,outputname,".csv"])

def csvmarkerremove(dataframeconcatenate,markerarray,outpfold, outputname):
    """
    Remove marker from phenograph input matrix
    :param dataframeconcatenate:
    :param markerarray:
    :return:
    """
    dataframeconcatenate = pd.read_csv(dataframeconcatenate,header=0)
    for i in range(len(markerarray)):
        dataframeconcatenate.drop(markerarray[i].strip('\n'), axis=1, inplace=True)
    dataframeconcatenate.to_csv("".join([outpfold, outputname, "_nomarker.csv"]), index=False, header=True)
    return "".join([outpfold, outputname, "_nomarker.csv"])

def fcs2csv(inputfcsfolder,outputcsvfolder):
    # set the extension value
    extension = 'fcs'
    # save path and input filename
    allfilenames = [i for i in glob.glob("".join([FCS_folder, '*.{}'.format(extension)]))]
    # save prefix name
    prefixname = []
    #
    for i in range(len(allfilenames)):
        sample = fk.Sample(allfilenames[i], subsample_count=None)
        sample.export_csv(source='raw', subsample=False,
                          filename="".join([allfilenames[i].split("/")[-1].split(".")[0],".csv"]),
                          directory=output_CSV_folder)

if __name__ == '__main__':
    print("Script start")
    raw_folder = pathchecker(options.input_folder)
    output_folder = pathchecker(options.output_folder)
    if options.analysis_name is None:
        print("Analysis name not provided. Please add -n flag.")
        sys.exit(1)
    else:
        pass
    if options.channel_excluded is None:
        removestep = False
    else:
        removestep = True
        marker = read_markers(options.channel_excluded)
    print(" ".join(["These markers will be removed:\n",'-'.join(marker)]))
    if options.f.upper() == "CSV":
        print("CSV files format are selected.")
        if removestep is False:
            print("remove step false")
        else:
            print("Concatenation of input CSV files")
            concfile = csvconcatenate(raw_folder, output_folder,options.analysis_name)
            print("Remove markers")
            pg_input = csvmarkerremove(concfile, marker,output_folder,options.analysis_name)
    elif options.f.upper() == "FCS":
        print("FCS files format are selected.")
        dirCSV = "".join([output_folder, "converted_fcs"])
        # create directory for the converted FCS to CSV
        createdir(dirCSV)
        # convert FCS to CSV
        fcs2csv(raw_folder,dirCSV)
        #
        if options.f.upper() == "CSV":
            print("CSV files format are selected.")
            if removestep is False:
                print("remove step false")
            else:
                print("Concatenation of input CSV files")
                concfile = csvconcatenate(raw_folder, output_folder, options.analysis_name)
                print("Remove markers")
                pg_input = csvmarkerremove(concfile, marker, output_folder, options.analysis_name)
    else:
        print("Please specify input format: csv or fcs")
        sys.exit(1)