import pandas as pd
import random
import optparse
import os
import glob

parser = optparse.OptionParser(usage=' -i [input csv folder] -o [output csv folder] -d [subsample lines]',
                               version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-d', action="store", dest="downsampleline", help='')
options, args = parser.parse_args()


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


# move to input folder
os.chdir(options.input_folder)
# create output folder
createdir(options.output_folder)
# read all csv file
all_files = glob.glob(os.path.join(options.input_folder, "*.csv"))
# get files name
names = [os.path.basename(x) for x in all_files]
# loop over all files and
for file_, name in zip(all_files, names):
    n = sum(1 for line in open(file_)) - 1
    s = int(options.downsampleline)
    skip = sorted(random.sample(range(1, n + 1), n - s))
    file_df = pd.read_csv(file_, header=0, skiprows=skip)
    file_df.to_csv("/".join([options.output_folder, "".join([name.split('.')[0], "_", options.downsampleline,
                                                             ".csv"])]), header=True, index=False)
