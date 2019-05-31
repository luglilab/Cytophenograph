import pandas as pd
import optparse
import os
#
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-f', action="store", dest="file_name", help='')
options, args = parser.parse_args()
#
def scaling(inputfilepath,filename,outputfolder):

    dfInput = pd.read_csv(inputfilepath, header=0, index_col="filename_cell")
    dfInput2= dfInput.copy()
    for i in list(dfInput2.columns.values):
        min_q = dfInput2[i].min()
        max_q = dfInput2[i].quantile(0.995)
        range_q = max_q - min_q
        dfInput2[i] = dfInput2[i].apply(lambda x: (x - min_q) / range_q)
    dfInput2.to_csv("".join([outputfolder, filename,"scaled.txt"]),header=True,index=True)
    print('Normalized 0 to 1 with MAX={}th percentile'.format(0.995*100),
          flush=True)

if __name__ == '__main__':
    scaling(inputfilepath="/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs/test_nomarker.csv",
            filename= os.path.splitext(os.path.basename("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs/test_nomarker.csv"))[0],
            outputfolder="/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/CD4_Phenograph_data/CD4_filtered_42_122/fcs/")

