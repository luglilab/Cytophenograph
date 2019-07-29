from PhenoFunctions import *

markerslist = "/home/spuccio/datadisk2/SP008_Phenograph_BMT/Jasper/Marker.txt"
outfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Jasper/test/30"
name = "test"
k = 30
Thread = 30
csvfolder = "/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/test"
infodf = filenameeventinfor(csvfolder)

pathmarkerfile, basenamemarkerfilepath = loadmarkers(markerslist)


markertoexclude = checkmarkers(markerfile="/".join([pathmarkerfile, basenamemarkerfilepath]),
                               data=infodf[2])




# dfPheno = runphenograph(marker=markertoexclude,
#                         concdf=infodf[1],
#                         kcoev=k,
#                         thread=Thread,
#                         outfold=outfolder,
#                         name=name)
#
# dfTsne = runtsne(marker=markertoexclude,
#                  concdf=infodf[1],
#                  thread=Thread,
#                  outfold=outfolder,
#                  name=name)
#
#
# dfAll = combine(data=infodf[2],
#                 tsne=dfTsne, pheno=dfPheno, outfold=outfolder, name=name)
#
# groupbycluster(alldf=dfAll, outfold=outfolder, name=name)
#
# groupbysample(alldf=dfAll, outfold=outfolder, name=name)
#
# validationplot(marker=markertoexclude, alldf=dfAll, outfold=outfolder, name=name)
#
# tsneplot(x=dfTsne.values, colors=dfPheno.values.flatten(),
#          outfold=outfolder, name=name)

