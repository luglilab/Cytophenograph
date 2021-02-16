from version import __version__
import optparse
from PhenoFunctions_v2_0 import *
import os
parser = optparse.OptionParser(usage='python ./Cytophenograph/cytophenograph.v2_0.py -i $abs_path/Cytophenograph/Test_dataset/CD8_Panel_II_channelvalues_GA_downSampled/ -o $abs_path/Cytophenograph/output_test -k 300 -m $abs_path/Cytophenograph/Test_dataset/CD8_bulk_markers_to_exclude.txt -n Test -t 10 -p $abs_path/Cytophenograph/Test_dataset/Info_file_bulk_Test.xlsx -c Parc', version='2.0')
parser.add_option('-i', action="store", dest="input_folder", help='Absolute path of folder with CSV files.')
parser.add_option('-o', action="store", dest="output_folder", help='Absolute path of output folder. TIPS: Please use an empty folder.')
parser.add_option('-k', action="store", dest="kmeancoef", help='Number for nearest neighbors search for Phenograph execution.')
parser.add_option('-m', action="store", dest="markerlist", help='Text file with features(channel name) to exclude during clustering execution.')
parser.add_option('-n', action="store", dest="analysis_name", help='Analysis name.')
parser.add_option('-t', action="store", dest="thread", help='Number of jobs.')
parser.add_option('-p', action="store", dest="pheno", help='Excel file with the following columns "Sample-Cell_type-EXP-ID-Time_point-Condition-Count", that will be integrated as metadata.')
parser.add_option('-c', type='choice', choices=['Phenograph', 'Parc', 'Both'],
                  dest="clustering",default = "Phenograph", help='Tool selecting for clustering. Both will execute Phenograph and Parc.')
parser.add_option("-s", action="store_true", dest="scale", default=False,help='Apply MinMaxScaler(1-100) before clustering.')
options, args = parser.parse_args()


if __name__ == '__main__':
    print("Script name: cytophenograph" )
    print("Script version:",__version__)
    print("Start")
    DictInfo = dict()

    run = Cytophenograph(info_file=options.pheno, input_folder=options.input_folder,
                         output_folder=options.output_folder,
                         k_coef=options.kmeancoef,
                         marker_list=options.markerlist,
                         analysis_name=options.analysis_name,
                         thread=int(options.thread),
                         tool= options.clustering,
                         scale=options.scale)
    DictInfo["Infofile"] = run.read_info_file()
    DictInfo["List_csv_files"] = run.import_all_event()
    DictInfo["adata_conc"]  = run.concatenate_dataframe(DictInfo["Infofile"],DictInfo["List_csv_files"])
    DictInfo["pathmarkerfile"] , DictInfo["basenamemarkerfilepath"]  = run.loadmarkers()
    DictInfo["markertoexclude"]  = run.checkmarkers(DictInfo["adata_conc"])
    if options.clustering == "Phenograph":
        print("Clustering tool selected is: Phenograph")
        DictInfo["phenograph_adata"] = run.runclustering(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["phenograph_adata"],"Phenograph")
        run.groupbysample(DictInfo["phenograph_adata"],"Phenograph")
        run.exporting(DictInfo["phenograph_adata"],options.clustering)
    elif options.clustering == "Parc":
        print("Clustering tool selected is: Parc")
        DictInfo["parc_adata"] = run.runparc(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["parc_adata"],"Parc")
        run.groupbysample(DictInfo["parc_adata"],"Parc")
        run.exporting(DictInfo["parc_adata"],options.clustering)
    elif options.clustering == "Both":
        print("Clustering tool selected is: Phenograph and Parc")
        DictInfo["phenograph_adata"] = run.runclustering(DictInfo["markertoexclude"], DictInfo["adata_conc"])
        run.groupbycluster(DictInfo["phenograph_adata"],"Phenograph")
        run.groupbysample(DictInfo["phenograph_adata"],"Phenograph")
        DictInfo["parc_adata"] = run.runparc(DictInfo["markertoexclude"], DictInfo["phenograph_adata"])
        run.groupbycluster(DictInfo["parc_adata"],"Parc")
        run.groupbysample(DictInfo["parc_adata"],"Parc")
        run.exporting(DictInfo["parc_adata"],options.clustering)
    print("End")
