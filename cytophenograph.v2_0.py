import optparse
from PhenoFunctions_v2_0 import *
import os
import pandas as pd
parser = optparse.OptionParser(usage='', version='1.0')
parser.add_option('-i', action="store", dest="input_folder", help='')
parser.add_option('-o', action="store", dest="output_folder", help='')
parser.add_option('-k', action="store", dest="kmeancoef", help='')
# parser.add_option('-k', action="append", dest='kmeancoef', default=[])
parser.add_option('-m', action="store", dest="markerlist", help='')
parser.add_option('-n', action="store", dest="analysis_name", help='')
parser.add_option('-t', action="store", dest="thread", help='')
parser.add_option('-p', action="store", dest="pheno", help='')
parser.add_option('-c', type='choice', choices=['Phenograph', 'Parc', 'Both'],
                  dest="clustering",default = "Phenograph", help='')
parser.add_option('-d', type='choice', choices=[True, False],dest="tsne",default = False, help='')
options, args = parser.parse_args()


if __name__ == '__main__':
    DictInfo = dict()

    run = Cytophenograph(info_file=options.pheno, input_folder=options.input_folder,
                         output_folder=options.output_folder,
                         k_coef=options.kmeancoef,
                         marker_list=options.markerlist,
                         analysis_name=options.analysis_name,
                         thread=int(options.thread),
                         tsne=options.tsne,
                         tool= options.clustering)
    #DictInfo["Log1"],DictInfo["Log2"] = run.create_logfile()
    #print(DictInfo["Log1"],DictInfo["Log2"])
    DictInfo["Infofile"] = run.read_info_file()
    DictInfo["List_csv_files"] = run.import_all_event()
    DictInfo["adata_conc"] = run.concatenate_dataframe(DictInfo["Infofile"],DictInfo["List_csv_files"])
    DictInfo["pathmarkerfile"] , DictInfo["basenamemarkerfilepath"]  = run.loadmarkers()
    DictInfo["markertoexclude"]  = run.checkmarkers(DictInfo["adata_conc"])
    if options.clustering == "Phenograph":
        print("Clustering tool selected is: Phenograph")
        DictInfo["phenograph_adata"] = run.runclustering(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["phenograph_adata"],"Phenograph")
        run.groupbysample(DictInfo["phenograph_adata"],"Phenograph")
        run.exporting(DictInfo["phenograph_adata"])
    elif options.clustering == "Parc":
        print("Clustering tool selected is: Parc")
        DictInfo["parc_adata"] = run.runparc(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["parc_adata"],"Parc")
        run.groupbysample(DictInfo["parc_adata"],"Parc")
        run.exporting(DictInfo["parc_adata"])
    elif options.clustering == "Both":
        print("Clustering tool selected is: Phenograph and Parc")
        DictInfo["phenograph_adata"] = run.runclustering(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["phenograph_adata"],"Phenograph")
        run.groupbysample(DictInfo["phenograph_adata"],"Phenograph")
        DictInfo["parc_adata"] = run.runparc(DictInfo["markertoexclude"], DictInfo["phenograph_adata"])
        run.groupbycluster(DictInfo["parc_adata"],"Parc")
        run.groupbysample(DictInfo["parc_adata"],"Parc")
        run.exporting(DictInfo["parc_adata"])
