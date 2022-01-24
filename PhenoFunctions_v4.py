import warnings
import anndata
import glob
import os
import sys
import pandas as pd
import phenograph as pg
import scanpy as sc
import parc
import umap
import logging
from flowsom import flowsom as flowsom
import tempfile
tmp = tempfile.NamedTemporaryFile()
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=300, facecolor='white',
                              figsize=(10, 10))
sc.settings.verbosity = 0
warnings.filterwarnings("ignore", category=FutureWarning)


class Cytophenograph:
    def __init__(self, info_file, input_folder, output_folder, k_coef, marker_list, analysis_name, thread,tool):
        self.info_file = info_file
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.k_coef = k_coef
        self.marker_list = marker_list
        self.analysis_name = analysis_name
        self.thread = thread
        self.tool = tool
        self.tmp_df = pd.DataFrame()
        self.adata = None
        self.adata_subset = None
        self.embedding = None
        self.palette = None
        self.marker = None
        self.markertoinclude = None
        self.marker_array = None
        self.anndata_list = []
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
        self.log.info("Part1: Files concatenation")
        # create empy list for save several df
        pandas_df_list = []
        # create list with anndata object

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
                                # ann_tmp.obs['Cell_type'] = cell_type.to_string().split(" ")[-1]
                                ann_tmp.obs['Cell_type'] = ''.join(e for e in cell_type.to_string() if e.isalnum())
                                #
                                exp = info_file['EXP'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                # ann_tmp.obs['EXP'] = exp.to_string().split(" ")[-1]
                                ann_tmp.obs['EXP'] = ''.join(e for e in exp.to_string() if e.isalnum())
                                #
                                id = info_file['ID'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                # ann_tmp.obs['ID'] = id.to_string().split(" ")[-1]
                                ann_tmp.obs['ID'] = ''.join(e for e in id.to_string() if e.isalnum())
                                #
                                time_point = info_file['Time_point'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Time_point'] = time_point.to_string().split(" ")[-1]
                                ann_tmp.obs['Time_point'] = ''.join(e for e in time_point.to_string() if e.isalnum())
                                #
                                condition = info_file['Condition'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                                # ann_tmp.obs['Condition'] = condition.to_string().split(" ")[-1]
                                ann_tmp.obs['Condition'] = ''.join(e for e in condition.to_string() if e.isalnum())
                                #
                                count = info_file['Count'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                                # ann_tmp.obs['Count'] = count.to_string().split(" ")[-1]
                                ann_tmp.obs['Count'] = ''.join(e for e in count.to_string() if e.isalnum())
                                self.anndata_list.append(ann_tmp)
                            else:
                                self.log.error("Error, this file {0} is not in the column Sample of Infofile. \n Please check sample name and Infofile".format(pandas_df_list[i].index[0][:-2]))
                                sys.exit(1)
                        tmp = self.anndata_list[0]
                        self.anndata_list.pop(0)
                        if len(self.anndata_list) ==1:
                            self.adata = tmp
                            self.adata.layers['raw_value'] = self.adata.X
                        else:
                            self.adata = tmp.concatenate(self.anndata_list)
                            self.adata.layers['raw_value'] = self.adata.X
                    except (ValueError, Exception):
                        self.log.error("Error. Please check Info File Header or CSV header.")
                        sys.exit(1)
        else:
            self.log.error("Error. Please check Info File Header or CSV header.")
            sys.exit(1)
        self.tmp_df = pd.DataFrame(self.adata.X, index=self.adata.obs.index)
        self.tmp_df.columns = self.adata.var_names
        return self.adata

    def loadmarkers(self):
        """
        Read marker filename with path
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

    def checkmarkers(self):
        """
        Check if marker in file is also a column of conc file
        :return:
        """
        # read marker file
        self.marker_array = [line.rstrip() for line in open(self.marker_list)]
        # read concatenate file
        for i in range(len(self.marker_array)):
            if self.marker_array[i] in self.adata.var_names.to_list():
                continue
            else:
                print("Marker {} not found in Matrix.".format(self.marker_array[i]))
                self.log.error("Marker {} not found in Matrix.".format(self.marker_array[i]))
                sys.exit(1)
        return self.marker_array

    def splitmarker(self):
        """

        """
        self.marker = self.adata.var_names.to_list()
        self.markertoinclude = [i for i in self.marker if i not in self.marker_array]
        return self.markertoinclude

    def runumap(self):
        """
        Function that perform UMAP dimensional reduction
        return: embedding
        """
        self.log.info("Part3: UMAP (Uniform Manifold Approximation and Projection) generation")
        reducer = umap.UMAP(random_state=42, n_neighbors=10, min_dist=0.5)
        embedding = reducer.fit_transform(self.adata_subset.X)
        return embedding

    def plot_umap(self):
        """
        Function per generation of pdf files with umap plot
        """
        self.createdir("/".join([self.output_folder, "".join(["Figures",self.tool])]))
        sc.settings.figdir = "/".join([self.output_folder, "".join(["Figures", self.tool])])
        if len(self.adata_subset.obs["pheno_leiden"].unique()) < 10:
            self.palette = "tab10"
        else:
            self.palette = "tab20"
        sc.pl.umap(self.adata_subset, color="pheno_leiden",
                   palette=self.palette, legend_fontoutline=2, show=False, add_outline=False, frameon=False,
                   title="UMAP Plot",
                   s=50, save=".".join(["".join([str(self.tool), "_cluster"]), "pdf"]))
        sc.pl.umap(self.adata_subset, color="pheno_leiden",
                   palette=self.palette, legend_fontoutline=2, show=False, add_outline=False, frameon=False,
                   legend_loc='on data', title="UMAP Plot",
                   s=50, save="_legend_on_data.".join(["".join([str(self.tool), "_cluster"]), "pdf"]))
        sc.pl.correlation_matrix(self.adata_subset, "pheno_leiden", show=False,
                                 save=".".join([self.tool, "pdf"]))
        for _ in list(self.adata_subset.var_names.unique()):
            sc.pl.umap(self.adata_subset, color=_, show=False, layer="scaled",
                       legend_fontoutline=1, na_in_legend=False, s=30,
                       title=_, palette='Viridis', groups=[_],
                       save=".".join([''.join(e for e in _ if e.isalnum()), "pdf"])
                       )

    def matrixplot(self):
        """
        Function for the generation of matrixplot sc.pl.matrixplot
        return:
        """
        sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                         dendrogram=True, vmin=-2, vmax=2, cmap='RdBu_r', layer="scaled",
                         show=False, swap_axes=False,return_fig=False,
                         save=".".join(["matrixplot_mean_z_score", "pdf"]))
        sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                         dendrogram=True, cmap='Blues', standard_scale='var',
                         colorbar_title='column scaled\nexpression', layer="scaled",
                         swap_axes=False,return_fig=False,
                         show=False,
                         save=".".join(["matrixplot_column_scaled_expression", "pdf"]))

    def runphenograph(self):
        """
        Function for execution of phenograph analysis
        :return:
        """
        self.log.info("Part2: Phenograph Clustering")
        self.log.info("Markers used for Phenograph clustering:")
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.log.info(self.adata_subset.var_names)
        self.log.info("Markers excluded for Phenograph clustering:")
        self.log.info(self.marker_array)
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                         zero_center=True, copy=True).X
        self.adata_subset.X = self.adata_subset.layers['scaled']
        self.adata_subset.obs['pheno_leiden'], graph, Q = pg.cluster(self.adata_subset.X, k=int(self.k_coef),
                                           seed=42,
                                           clustering_algo="leiden",
                                           directed=False,
                                           prune=False, min_cluster_size=1,
                                           n_jobs=int(self.thread))
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden']
        self.adata.obs['Phenograph_cluster'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.embedding = self.runumap()
        self.adata.obsm['X_umap'] = self.embedding
        self.adata_subset.obsm['X_umap'] = self.embedding
        self.tmp_df = pd.DataFrame(self.adata.X, columns=self.adata.var_names)
        self.tmp_df['UMAP_1'] = self.embedding[:, 0]
        self.tmp_df['UMAP_2'] = self.embedding[:, 1]
        self.tmp_df['Cluster_Phenograph'] = self.adata_subset.obs['pheno_leiden']
        self.plot_umap()
        self.matrixplot()
        self.tmp_df.to_csv("/".join([self.output_folder, "_ConcatenatedCells.".join(["_".join([self.analysis_name]), "csv"])]),
                           header=True, index=False)
        return self.adata

    def runparc(self):
        """
        function for execution of
        :return:
        """
        self.log.info("Part2: PARC Clustering")
        self.log.info("Markers used for Parc clustering:")
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.log.info(self.adata_subset.var_names)
        self.log.info("Markers excluded for Phenograph clustering:")
        self.log.info(self.marker_array)
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                         zero_center=True, copy=True).X
        self.adata_subset.X = self.adata_subset.layers['scaled']
        p = parc.PARC(self.adata_subset.X, random_seed=42, knn=int(self.k_coef),
                      jac_std_global='median', jac_weighted_edges=False,
                      small_pop=1,
                      num_threads=int(self.thread))
        p.run_PARC()
        self.adata_subset.obs['pheno_leiden'] = [str(i) for i in p.labels]
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden']
        self.adata.obs['Parc_cluster'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.embedding = self.runumap()
        self.adata.obsm['X_umap'] = self.embedding
        self.adata_subset.obsm['X_umap'] = self.embedding
        self.tmp_df = pd.DataFrame(self.adata.X, columns=self.adata.var_names)
        self.tmp_df['UMAP_1'] = self.embedding[:, 0]
        self.tmp_df['UMAP_2'] = self.embedding[:, 1]
        self.tmp_df['Cluster_Parc'] = self.adata_subset.obs['pheno_leiden']
        self.plot_umap()
        self.matrixplot()
        self.tmp_df.to_csv("/".join([self.output_folder, "_ConcatenatedCells.".join(["_".join([self.analysis_name]), "csv"])]),
                           header=True, index=False)
        return self.adata

    def runflowsom(self):
        """
        function for execution of
        :return:
        """
        self.log.info("Part2: Flowsom Clustering")
        self.log.info("Markers used for Flowsom clustering:")
        self.log.info(self.markertoinclude)
        self.log.info("Markers excluded for Flowsom clustering:")
        self.log.info(self.marker_array)
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                         zero_center=True, copy=True).X
        self.adata.to_df().to_csv(tmp.name, header=True, index=False)
        tt = flowsom(tmp.name, if_fcs=False,if_drop=True,
                     drop_col=self.markertoinclude)
        sample_df = tt.df
        tt.som_mapping(50, 50, tt.df.shape[1],
                       sigma =2.5,
                       lr=0.1,
                       batch_size=100,
                       neighborhood='gaussian',
                       if_fcs=False,
                       seed=10)
        from sklearn.cluster import AgglomerativeClustering
        tt.meta_clustering(AgglomerativeClustering,
                           # cluster_class: e.g. KMeans, a cluster class, like "from sklearn.cluster import KMeans"
                           5,  # min_n: e.g. 10, the min proposed number of clusters
                           31,  # max_n: e.g. 31, the max proposed number of clusters
                           10,  # iter_n: e.g 10, the iteration times for each number of clusters
                           resample_proportion=0.6,
                           # resample_proportion: e.g 0.6, the proportion of re-sampling when computing clustering
                           verbose=False  # verbose: e.g. False, whether print out the clustering process
                           )
        tt.labeling()
        # output_df = tt.df  # new column added: category
        output_tf_df = tt.tf_df  # new column added: category
        output_tf_df.set_index(self.adata.obs.index, inplace=True)
        TMP = output_tf_df['category'].astype(int).sort_values().unique()
        res_dct = {TMP[i]: i for i in range(len(TMP))}
        output_tf_df['category'] = output_tf_df['category'].map(res_dct)
        self.adata_subset.obs['pheno_leiden'] = output_tf_df['category'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden']
        self.adata.obs['Cluster_Flowsom'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata_subset.X = self.adata_subset.layers['scaled']
        self.embedding = self.runumap()
        self.adata.obsm['X_umap'] = self.embedding
        self.adata_subset.obsm['X_umap'] = self.embedding
        self.tmp_df = pd.DataFrame(self.adata.X, columns=self.adata.var_names)
        self.tmp_df['UMAP_1'] = self.embedding[:, 0]
        self.tmp_df['UMAP_2'] = self.embedding[:, 1]
        self.tmp_df['Cluster_Flowsom'] = self.adata_subset.obs['pheno_leiden']
        self.plot_umap()
        self.matrixplot()
        self.tmp_df.to_csv(
            "/".join([self.output_folder, "_ConcatenatedCells.".join(["_".join([self.analysis_name]), "csv"])]),
            header=True, index=False)
        return self.adata

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

    def groupbycluster(self):
        """
        Function for generation of csv with different clusters
        :parm tool:
        :return:
        """
        # make dir
        self.createdir("/".join([self.output_folder, "".join(["FCScluster", self.tool])]))
        self.adata.obs['cluster'] = self.adata.obs['cluster'].astype(int)
        for _ in range(self.adata.obs['cluster'].unique().min(), self.adata.obs['cluster'].unique().max() + 1):
            tmp = self.adata[self.adata.obs['cluster'].isin([_])].to_df()
            tmp['UMAP_1'] = self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:,0]
            tmp['UMAP_2'] = self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:,1]
            if (self.tool == "Phenograph"):
                tmp['Phenograph'] = _
            elif (self.tool == "Parc"):
                tmp['Parc'] = _
            else:
                tmp['Flowsom'] = _-1
            tmp.to_csv("/".join([self.output_folder, "".join(["FCScluster",
                                                              self.tool]), "".join([self.analysis_name, "_",
                                                                               str(_), ".csv"])]),
                       header=True, index=False)

    def groupbysample(self):
        """
        Function for generation of csv with different clusters
        """
        self.adata.obs['cluster'] = self.adata.obs['cluster'].astype(int)
        # make dir
        self.createdir("/".join([self.output_folder,"".join(["FCSsample",self.tool])  ]))
        self.createdir("/".join([self.output_folder, "".join(["Cluster_frequency", self.tool])]))
        # tmp df
        tmp = self.adata.to_df()
        # change index
        tmp.set_index(self.adata.obs['Sample'], inplace=True)
        # create columns
        tmp['UMAP_1'] = self.adata.obsm['X_umap'][:,0]
        tmp['UMAP_2'] = self.adata.obsm['X_umap'][:,1]
        tmp[self.tool] = self.adata.obs['cluster'].values
        tmp["cluster"] = self.adata.obs['cluster'].values
        # get unique filenames
        unique_filename = self.adata.obs['Sample'].unique()
        # get unique number of cluster
        unique_Phenograph = self.adata.obs['cluster'].unique()
        #
        dfCounts = pd.DataFrame(index=range(min(unique_Phenograph), max(unique_Phenograph) + 1))
        # generate Tot_percentage file
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].cluster.value_counts(
                normalize=True).reindex(tmp.cluster.unique(), fill_value=0)
        # compute percentage
        dfCounts = dfCounts * 100
        # save
        dfCounts.index.name = 'Cluster'
        dfCounts.to_csv("/".join(["/".join([self.output_folder, "".join(["Cluster_frequency", self.tool])]),
                                  "".join(["Tot_percentage", ".csv"])]))
        # create empty dataframe
        dfCounts = pd.DataFrame(index=range(min(unique_Phenograph), max(unique_Phenograph) + 1))
        # generate Tot_counts file
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].cluster.value_counts().reindex(
                tmp.cluster.unique(), fill_value=0)
        # save
        dfCounts.index.name = 'Cluster'
        dfCounts.to_csv("/".join(["/".join([self.output_folder, "".join(["Cluster_frequency", self.tool])]), "".join(["Tot_counts", ".csv"])]))
        del tmp['cluster']
        # save samples
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].to_csv(
                "/".join([self.output_folder, "".join(["FCSsample", self.tool]), "".join([str(unique_filename[i]), "_", self.analysis_name,
                                                                    ".csv"])]),
                header=True, index=False)

    def exporting(self):
        """
        Export to h5ad file.
        """
        self.log.info("Part4: Output Generation")
        old_names = self.adata.var_names
        new_names = []
        for _ in range(len(old_names)):
            if old_names[_].startswith("Comp-"):
                new_names.append(old_names[_].split(":: ")[-1])
            else:
                new_names.append(old_names[_])
        self.adata.var = pd.DataFrame(old_names, new_names)
        del self.adata.var[0]
        self.adata.var['original_names'] = old_names
        if self.tool == "Phenograph" or self.tool == "Parc":
            self.adata.obs[self.tool+"_"+str(self.k_coef)] = self.adata.obs['cluster'].astype("str")
            del self.adata.obs['cluster']
            del self.adata.obs[self.tool+"_"+str(self.k_coef)]
            self.adata.obs["".join([str(self.tool), "_cluster"])] = self.adata.obs["".join([str(self.tool), "_cluster"])].astype('category')
            self.adata.X = sc.pp.scale(self.adata, max_value=6, zero_center=True, copy=True).X
            self.adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name, self.k_coef]), "h5ad"])]))
        elif self.tool == "Flowsom":
            del self.adata.obs['cluster']
            self.adata.X = sc.pp.scale(self.adata, max_value=6,zero_center=True, copy=True).X
            self.adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "h5ad"])]))

