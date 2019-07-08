import optparse
from PhenoFunctions import *
import os
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-k', action="store", dest="kmeancoef", help='')
parser.add_option('-m', action="store", dest="markerlist", help='')
parser.add_option('-n', action="store", dest="analysis_name", help='')
parser.add_option('-t', action="store", dest="thread", help='')
#parser.add_option('-f', type='choice', choices=['csv', 'fcs', 'FCS', 'CSV'], help='')
options, args = parser.parse_args()


markerslist = "/home/spuccio/datadisk2/SP008_Phenograph_BMT/Jasper/Marker.txt"
outfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Jasper/test/30"
name = "test"
k = 30
Thread = 30
#csvfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/test"


if __name__ == '__main__':
    infodf = filenameeventinfor(options.input_folder)

    pathmarkerfile, basenamemarkerfilepath = loadmarkers(options.markerlist)

    markertoexclude = checkmarkers(markerfile="/".join([pathmarkerfile, basenamemarkerfilepath]),
                                   data=infodf[2])


    dfPheno = runphenograph(marker=markertoexclude,
                            concdf=infodf[1],
                            kcoev=int(options.kmeancoef),
                            thread=int(options.thread),
                            outfold=options.output_folder,
                            name=options.analysis_name)

    dfTsne = runtsne(marker=markertoexclude,
                     concdf=infodf[1],
                     thread=int(options.thread),
                     outfold=options.output_folder,
                     name=options.analysis_name)

    dfAll = combine(data=infodf[2],
                    tsne=dfTsne, pheno=dfPheno, outfold=options.output_folder, name=options.analysis_name)

    groupbycluster(alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)

    groupbysample(alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)

    validationplot(marker=markertoexclude, alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)
    x = dfTsne[['Tsne_1','Tsne_2']]
    tsneplot(x=dfTsne.values, colors=dfPheno["Phenograph"], outfold=options.output_folder, name=options.analysis_name)
