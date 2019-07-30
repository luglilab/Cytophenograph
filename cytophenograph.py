import optparse
from PhenoFunctions import *
import os
import pandas as pd
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-k', action="store", dest="kmeancoef", help='')
#parser.add_option('-k', action="append", dest='kmeancoef', default=[])
parser.add_option('-m', action="store", dest="markerlist", help='')
parser.add_option('-n', action="store", dest="analysis_name", help='')
parser.add_option('-t', action="store", dest="thread", help='')
#parser.add_option('-f', type='choice', choices=['csv', 'fcs', 'FCS', 'CSV'], help='')
options, args = parser.parse_args()


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

    dfUmap = runumap(marker=markertoexclude,
                     concdf=infodf[1],
                     thread=int(options.thread),
                     outfold=options.output_folder,
                     name=options.analysis_name)

    dfAll = combine(data=infodf[2],
                    tsne=dfTsne, pheno=dfPheno, outfold=options.output_folder, name=options.analysis_name)

    dfAll.to_csv("/".join([options.output_folder, ".".join([options.analysis_name,"_concatenad.txt"])]), sep="\t",
                 header=True, index=True)

    groupbycluster(alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)

    groupbysample(alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)

    validationplot(marker=markertoexclude, alldf=dfAll, outfold=options.output_folder, name=options.analysis_name)

    tsne_umap_plot(x=dfTsne, colors=dfPheno, outfold=options.output_folder, name=options.analysis_name,
             kind="Tsne")
    tsne_umap_plot(x=dfUmap, colors=dfPheno, outfold=options.output_folder, name=options.analysis_name,
             kind="Umap")
