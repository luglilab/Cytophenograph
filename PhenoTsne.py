from PhenoFunctions import *

markerslist = "/home/spuccio/datadisk2/SP008_Phenograph_BMT/Jasper/Marker.txt"
outfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Jasper/CD4_Ag_Phenograph//30"
name = "CD4_Ag_allfiles_k30"
k = 30
Thread = 30
csvfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Exported_FJ10_CD4_Ag"
infodf = filenameeventinfor(csvfolder)

pathmarkerfile, basenamemarkerfilepath = loadmarkers(markerslist)


markertoexclude = checkmarkers(markerfile="/".join([pathmarkerfile, basenamemarkerfilepath]),
                               data=infodf[2])

dfPheno = runphenograph(marker=markertoexclude,
                        concdf=infodf[1],
                        kcoev=k,
                        thread=Thread,
                        outfold=outfolder,
                        name=name)

dfTsne = runtsne(marker=markertoexclude,
                 concdf=infodf[1],
                 thread=Thread,
                 outfold=outfolder,
                 name=name)


dfAll = combine(data=infodf[2],
                tsne=dfTsne, pheno=dfPheno, outfold=outfolder, name=name)

groupbycluster(alldf=dfAll, outfold=outfolder, name=name)

groupbysample(alldf=dfAll, outfold=outfolder, name=name)


