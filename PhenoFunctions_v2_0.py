from version import __version__
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import anndata
import glob
import os
import sys
import numpy as np
import pandas as pd
import re
import phenograph as pg
import scanpy as sc
import parc
import umap
import logging
from sklearn import preprocessing
import matplotlib.pyplot as plt


class Cytophenograph:
    def __init__(self, info_file, input_folder, output_folder, k_coef, marker_list, analysis_name, thread,tool,scale):
        self.info_file = info_file
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.k_coef = k_coef
        self.marker_list = marker_list
        self.analysis_name = analysis_name
        self.thread = thread
        self.tool = tool
        self.tmp_df = pd.DataFrame()
        self.scale = scale
        # sys.stdout = open("/".join([self.output_folder,"log.txt"]), 'a')
        self.log = logging.getLogger()
        self.log.setLevel(logging.INFO)
        format = logging.Formatter("%(asctime)s %(threadName)-11s %(levelname)-10s %(message)s")
        #
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        self.log.addHandler(ch)
        # 
        fh = logging.FileHandler("/".join([self.output_folder,"log.txt"]), "w")
        fh.setFormatter(format)
        self.log.addHandler(fh)
        self.log.info("Name of this analysis: {}".format(marker_list))
        self.log.info("Input folder: {}".format(input_folder))
        self.log.info("Output folder: {}".format(output_folder))
        self.log.info("Info file: {}".format(info_file))
        self.log.info("Phenograph K-coef: {}".format(k_coef))
        self.log.info("Marker list file: {}".format(marker_list))
        self.log.info("Clustering tool option: {}".format(tool))
        
    def read_info_file(self):
        """
        Read info file methods
        :return: pandas dataframe
        """
        df_info = pd.read_excel(self.info_file, header=0)
        return df_info

    def import_all_event(self):
        """
        scan csv folder, save csv files names in list
        :return: list
        """
        # change folder
        os.chdir(self.input_folder)
        # create array with all csvfile name and path
        all_files = glob.glob(os.path.join(self.input_folder, "*.csv"))
        # get files name
        names = [os.path.basename(x) for x in all_files]
        #
        list_with_file_name_and_path = []
        for file_, name in zip(all_files, names):
            list_with_file_name_and_path.append(file_)
        return list_with_file_name_and_path

    def create_df(self, path_csv_file):
        """
        create dataframe file csv
        :return:
        """
        df = pd.read_csv(path_csv_file, header=0)
        barcode = []
        names = os.path.basename(path_csv_file)
        num_lines = df.shape[0]
        for _ in range(1, num_lines + 1):
            barcode.append("_".join([names.split(".")[0], str(_)]))
        df.index = barcode
        return df

    def concatenate_dataframe(self,info_file, csv_list):
        """

        :param csv_list:
        :return:
        """
        # create empy list for save several df
        pandas_df_list = []
        # create list with anndata object
        anndata_list = []
        # loop over csv file name
        for i in range(len(csv_list)):
            # append df in pandas_df_list
            pandas_df_list.append(self.create_df(csv_list[i]))
        # check header
        if all([len(pandas_df_list[0].columns.intersection(df.columns)) == pandas_df_list[0].shape[1]
                for df in pandas_df_list]):
                    try:
                        for i in range(len(pandas_df_list)):
                            # print(pandas_df_list[i].index[0][:-2]) sample id derivated
                            # save column with Sample name in list
                            Sample_list = info_file["Sample"].tolist()
                            # check if Sample name are in the anndata index
                            if pandas_df_list[i].index[0][:-2] in Sample_list:
                                ann_tmp = anndata.AnnData(pandas_df_list[i])
                                ann_tmp.obs['Sample'] = pandas_df_list[i].index[0][:-2]
                                #
                                cell_type = info_file['Cell_type'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Cell_type'] = cell_type.to_string().split(" ")[-1]
                                #
                                exp = info_file['EXP'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['EXP'] = exp.to_string().split(" ")[-1]
                                #
                                id = info_file['ID'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['ID'] = id.to_string().split(" ")[-1]
                                #
                                time_point = info_file['Time_point'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Time_point'] = time_point.to_string().split(" ")[-1]
                                #
                                condition = info_file['Condition'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Condition'] = condition.to_string().split(" ")[-1]
                                #
                                count = info_file['Count'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Count'] = count.to_string().split(" ")[-1]
                                anndata_list.append(ann_tmp)
                            else:
                                self.log.error("Error, this file {0} is not in the column Sample of Infofile. \n Please check sample name and Infofile".format(pandas_df_list[i].index[0][:-2]))
                                sys.exit(1)
                        tmp = anndata_list[0]
                        anndata_list.pop(0)
                        adata_conc = tmp.concatenate(anndata_list)
                    except (ValueError, Exception):
                        self.log.error("Error. Please check Info File Header or CSV header.")
                        sys.exit(1)
        self.tmp_df = pd.DataFrame(adata_conc.X,index=adata_conc.obs.index)
        self.tmp_df.columns = adata_conc.var_names
        pd.merge(self.tmp_df,adata_conc.obs,left_index=True, right_index=True).to_csv("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "txt"])]),header=True, index=False)
        return adata_conc

    def loadmarkers(self):
        """
        Read marker filename with path
        :param markerfile:
        :return: array with filepath and filename
        """
        if os.path.exists(self.marker_list) is True:
            markerfilepath = os.path.split(self.marker_list)[0]
            markerfilename = os.path.split(self.marker_list)[1]
            return markerfilepath, markerfilename
        else:
            print("File does not exist. Exit")
            self.log.error("File does not exist. Exit")
            sys.exit(1)

    def checkmarkers(self, data):
        """
        Check if marker in file is also a column of conc file
        :param markerfile:
        :param data:
        :return:
        """
        # read marker file
        marker_array = [line.rstrip() for line in open(self.marker_list)]
        # read concatenate file
        for i in range(len(marker_array)):
            if marker_array[i] in data.var_names.to_list():
                continue
            else:
                print("Marker {} not found in Matrix.".format(marker_array[i]))
                self.log.error("Marker {} not found in Matrix.".format(marker_array[i]))
                sys.exit(1)
        return marker_array

    def runclustering(self,markertoexclude, adata):
        """
        Function for execution of phenograph analysis
        :param markertoexclude:
        :param adata:
        :return:
        """
        marker = adata.var_names.to_list()
        markertoinclude = [i for i in marker if i not in markertoexclude]
        data = adata[:, markertoinclude].to_df()
        self.log.info("Markers used for Phenograph clustering:")
        self.log.info(data.columns)
        if (self.scale == True):
            min_max_scaler = preprocessing.MinMaxScaler((1, 100))
            x_scaled = min_max_scaler.fit_transform(data.values)
            data = pd.DataFrame(x_scaled,
            columns=data.columns)
        self.new_head=[]
        self.new_head.append([column.split("::")[-1] for column in data])
        data.columns = self.new_head
        data.diff().hist(color="k", alpha=0.5, bins=50, grid=False, xlabelsize=8,ylabelsize=8)
        plt.tight_layout()
        plt.savefig("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "pdf"])]))
        communities, graph, Q = pg.cluster(data.values, k=int(self.k_coef),directed=False,
        prune=False,min_cluster_size=1,n_jobs=int(self.thread))
        # create dataframe with Phenograph output
        self.dfPheno = pd.DataFrame(communities)
        # shift of one unit the name of cluster
        self.dfPheno["Phenograph"] = self.dfPheno[0] + 1
        # remove first column
        self.dfPheno = self.dfPheno.drop(columns=[0], axis=1)
        self.dfPheno.set_index(adata.obs.index, inplace=True)
        adata.obs['cluster'] = self.dfPheno
        adata.obs['Phenograph_cluster'] = self.dfPheno
        reducer = umap.UMAP(random_state=42, n_neighbors=10, min_dist=0.001)
        embedding = reducer.fit_transform(data.values)
        adata.obsm['X_umap'] = embedding
        self.tmp_df = self.tmp_df.astype(int)
        self.tmp_df['UMAP_1'] = embedding[:,0]
        self.tmp_df['UMAP_2'] = embedding[:,1]
        self.tmp_df['Cluster_Phenograph'] = self.dfPheno
        self.tmp_df.to_csv("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "csv"])]),header=True, index=False)
        return adata

    def runparc(self, markertoexclude, adata):
        """
        function for execution of
        :param adata:
        :param kcoev:
        :param thread:
        :return:
        """
        marker = adata.var_names.to_list()
        markertoinclude = [i for i in marker if i not in markertoexclude]
        data = adata[:, markertoinclude].to_df()
        self.log.info("Markers used for PARC clustering:")
        self.log.info(data.columns)
        if (self.scale == True):
            min_max_scaler = preprocessing.MinMaxScaler((1, 100))
            x_scaled = min_max_scaler.fit_transform(data.values)
            data = pd.DataFrame(x_scaled,
            columns=data.columns)
        self.new_head=[]
        self.new_head.append([column.split("::")[-1] for column in data])
        data.columns = self.new_head
        data.diff().hist(color="k", alpha=0.5, bins=50, grid=False, xlabelsize=8,ylabelsize=8)
        plt.tight_layout()
        plt.savefig("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "pdf"])]))
        p = parc.PARC(data.values, random_seed=42,
                      jac_std_global='median',
                      small_pop = 100,
                      num_threads=int(self.thread))
        p.run_PARC()
        adata.obs['cluster'] = [str(i) for i in p.labels]
        adata.obs['Parc_cluster'] = [str(i) for i in p.labels]
        reducer=umap.UMAP(random_state=42,
                            n_neighbors=10,
                            min_dist=0.001)
        embedding = reducer.fit_transform(data.values)
        adata.obsm['X_umap'] = embedding
        self.tmp_df = self.tmp_df.astype(int)
        self.tmp_df['UMAP_1'] = embedding[:,0]
        self.tmp_df['UMAP_2'] = embedding[:,1]
        if (self.tool == "Both"):
            self.tmp_df['Cluster_Phenograph'] = self.dfPheno
        self.tmp_df['Cluster_PARC'] = [str(i) for i in p.labels]
        self.tmp_df.to_csv("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "csv"])]),header=True, index=False)
        return adata

    def createdir(self,dirpath):
        """
        Make dir function and check if directory is already exists
        :param dirpath: string with path and directory name
        :return:
        """
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)
            print(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
            self.log.info(" ".join(["Directory", dirpath.split("/")[-1], "Created"]))
        else:
            print(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))
            self.log.info(" ".join(["Directory", dirpath.split("/")[-1], "already exists"]))

    def groupbycluster(self, adata,tool):
        """
        Function for generation of csv with different clusters
        :param adata:
        :return:
        """
        # make dir
        self.createdir("/".join([self.output_folder, "".join(["FCScluster",tool])]))
        adata.obs['cluster'] = adata.obs['cluster'].astype(int)
        for _ in range(adata.obs['cluster'].unique().min(), adata.obs['cluster'].unique().max() + 1):
            tmp = adata[adata.obs['cluster'].isin([_])].to_df()
            tmp['UMAP_1'] = adata[adata.obs['cluster'].isin([_])].obsm['X_umap'][:,0]
            tmp['UMAP_2'] = adata[adata.obs['cluster'].isin([_])].obsm['X_umap'][:,1]
            if (self.tool == "Phenograph"):
                tmp['Phenograph'] = _
            elif (self.tool == "Parc"):
                tmp['Parc'] = _
            else:
                tmp['Cluster'] = _
            tmp.to_csv("/".join([self.output_folder, "".join(["FCScluster",tool]), "".join([self.analysis_name, "_", str(_), ".csv"])]),
                       header=True, index=False)

    def groupbysample(self, adata, tool):
        """
        Function for generation of csv with different clusters
        :param alldf: Combine output
        file_name,event,FSC-A,FSC-H,FSC-W,SSC-A,SSC-H,Comp-APC-A,......,Tbet,Time,,Tsne_1,Tsne_2,Phenograph
        :param outfold: output folder
        :param name: output name analysis
        :return: None
        """
        adata.obs['cluster'] = adata.obs['cluster'].astype(int)
        # make dir
        self.createdir("/".join([self.output_folder,"".join(["FCSsample",tool])  ]))
        # tmp df
        tmp = adata.to_df()
        # change index
        tmp.set_index(adata.obs['Sample'], inplace=True)
        # create columns
        tmp['UMAP_1'] = adata.obsm['X_umap'][:,0]
        tmp['UMAP_2'] = adata.obsm['X_umap'][:,1]
        if (self.tool == "Phenograph"):
            tmp["Phenograph"] = adata.obs['cluster'].values
        elif (self.tool == "Parc"):
            tmp["Parc"] = adata.obs['cluster'].values
        else:
            tmp["cluster"] = adata.obs['cluster'].values
        tmp["cluster"] = adata.obs['cluster'].values
        # get unique filenames
        unique_filename = adata.obs['Sample'].unique()
        # get unique number of cluster
        unique_Phenograph = adata.obs['cluster'].unique()
        #
        dfCounts = pd.DataFrame(index=range(min(unique_Phenograph), max(unique_Phenograph) + 1))
        # generate Tot_percentage file
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].cluster.value_counts(
                normalize=True).reindex(tmp.cluster.unique(), fill_value=0)
        # compute percentage
        dfCounts = dfCounts * 100
        # save
        dfCounts.to_csv("/".join([self.output_folder, "".join(["FCSsample",tool]), "".join(["Tot_percentage", ".txt"])]))
        # create empty dataframe
        dfCounts = pd.DataFrame(index=range(min(unique_Phenograph), max(unique_Phenograph) + 1))
        # generate Tot_counts file
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].cluster.value_counts().reindex(
                tmp.cluster.unique(), fill_value=0)
        # save
        dfCounts.to_csv("/".join([self.output_folder, "".join(["FCSsample",tool]), "".join(["Tot_counts", ".txt"])]))
        del tmp['cluster']
        # save samples
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].to_csv(
                "/".join([self.output_folder, "".join(["FCSsample",tool]), "".join([str(unique_filename[i]), "_", self.analysis_name,
                                                                    ".csv"])]),
                header=True, index=False)

    def exporting(self, adata,tool):
        """
        Export to h5ad file. 
        """
        old_names = adata.var_names
        new_names = []
        for _ in range(len(old_names)):
            if old_names[_].startswith("Comp-"):
                new_names.append(old_names[_].split(":: ")[-1])
            else:
                new_names.append(old_names[_])
        adata.var = pd.DataFrame(old_names, new_names)
        del adata.var[0]
        adata.var['original_names'] = old_names
        if (tool == "Phenograph"):
            adata.obs[tool+"_"+str(self.k_coef)] = adata.obs['cluster'].astype("str")
            del adata.obs['cluster']
            del adata.obs[tool+"_"+str(self.k_coef)]
            adata.obs['Phenograph_cluster'] = adata.obs['Phenograph_cluster'].astype('category')
            adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name, self.k_coef]), "h5ad"])]))
        elif (tool == "Parc"):
            adata.obs['Parc_cluster'] = adata.obs['cluster'].astype("str")
            adata.obs['Parc_cluster'] = adata.obs['Parc_cluster'].astype('category') 
            del adata.obs['cluster']
            adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "h5ad"])]))
        elif (tool == "Both"):
            adata.obs["Phenograph_cluster"] = self.dfPheno.astype("str")
            del adata.obs['cluster']
            adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "h5ad"])]))
            
