import warnings
import anndata
import glob
import os
import sys
import pandas as pd
import phenograph as pg
import scanpy as sc
import pyVIA.core as via
import umap
import logging
import pickle
import tempfile
import matplotlib
import scanorama
import scipy
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import fcsy
from fcsy import DataFrame
import subprocess
from scipy.cluster.hierarchy import dendrogram, linkage
matplotlib.use('Agg')
import seaborn as sns
warnings.filterwarnings('ignore')
# import traceback
import numpy as np
tmp = tempfile.NamedTemporaryFile()
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi = 300, facecolor = 'white', dpi_save = 330,
                              figsize = (10, 10))
sc.settings.verbosity = 0
warnings.filterwarnings("ignore", category = FutureWarning)
from palette import palette28,palette102
import scprep
class CustomFormatter(logging.Formatter):
    FORMATS = {
        logging.INFO: "###%(msg)s",
        logging.WARNING: "$$$%(msg)s",
        logging.ERROR: "@@@%(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


class Cytophenograph:
    def __init__(self, info_file, input_folder, output_folder, k_coef, marker_list, analysis_name, thread, tool, batch,
                 batchcov, mindist, spread, runtime, knn, resolution, maxclus, downsampling, cellnumber,
                 filetype,arcsinh):
        self.info_file = info_file
        self.input_folder = input_folder
        self.output_folder = output_folder
        # self.k_coef = k_coef
        self.marker_list = marker_list
        # analysis file name 
        import re
        self.analysis_name = re.sub('[\W]+', '', analysis_name.replace(' ', '_'))
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
        self.flowsomshape = None
        self.anndata_list = []
        self.outfig = None
        self.tmp = None
        self.dpi = 300
        self.fileformat = "pdf"  # insert svg to change figure format
        self.newheader = []
        self.n_neighbors = 5
        self.log = logging.getLogger()
        self.log.setLevel(logging.INFO)
        self.scanorama = batch
        self.batchcov = batchcov
        self.runtime = runtime
        self.cleaning = {}
        self.target_cells = 0.1
        if self.tool == "Phenograph":
            self.k_coef = k_coef
        if self.tool == "VIA":
            self.knn = knn
            self.resolution = resolution
        if self.tool == "FlowSOM":
            self.maxclus = str(maxclus)
            self.flowsomDF = pd.DataFrame()
        self.listmarkerplot = None
        self.concatenate_fcs = None
        self.path_flowai = os.path.dirname(os.path.realpath(__file__)) + '/flowai.Rscript'
        self.path_flowsom = os.path.dirname(os.path.realpath(__file__)) + '/flowsom.Rscript'
        self.mindist = float(mindist)
        self.spread = float(spread)
        self.downsampling = downsampling
        self.cellnumber = cellnumber
        self.filetype = filetype
        if self.filetype == "FCS":
            self.arcsinh = arcsinh
        else:
            self.arcsinh = False
        self.root_user = [1]
        self.fnull = open(os.devnull, 'w')
        ch = logging.StreamHandler()
        ch.setFormatter(CustomFormatter())
        self.log.addHandler(ch)
        self.palette28 = palette28
        self.palette102 = palette102

        if self.runtime == 'UMAP': self.tool = 'UMAP'
        self.log.info("Runtime: {}".format(self.runtime))
        self.log.info("DownSampling: {}".format(self.downsampling))
        if self.downsampling != 'All':
            self.log.info(" >> Event Number: {}".format(self.cellnumber))
        self.log.warning("PART 1")

    def read_fcs(self, path_csv_file):
        """
        Read FCS file version 3 and convert in pandas dataframe
        Returns: Pandas Dataframe
        """
        df = DataFrame.from_fcs(path_csv_file, channel_type = 'multi')
        df.columns = df.columns.map(' :: '.join)
        df.columns = df.columns.str.replace('[\",\']', '')
        df.columns = df.columns.str.replace(' $', '', regex = True)
        df.columns = df.columns.str.replace(" \:: $", "", regex = True)
        if self.downsampling == "Balanced":
            if self.cellnumber < df.shape[0]:
                df = df.sample(n = int(self.cellnumber), random_state = 42,
                               ignore_index = True)
            else:
                pass
        barcode = []
        names = os.path.basename(path_csv_file)
        num_lines = df.shape[0]
        for _ in range(1, num_lines + 1):
            barcode.append("_".join([names.split(".")[0], str(_)]))
        df.index = barcode
        return df

    def read_info_file(self):
        """
        Read info file methods
        :return: pandas dataframe
        """
        df_info = pd.read_excel(self.info_file, header = 0)
        return df_info

    def import_all_event(self):
        """
        scan csv folder, save csv files names in list
        :return: list
        """
        # change folder
        os.chdir(self.input_folder)
        # create array with all csvfile name and path
        if self.filetype == "CSV":
            all_files = glob.glob(os.path.join(self.input_folder, "*.csv"))
        else:
            all_files = glob.glob(os.path.join(self.input_folder, "*.fcs"))
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
        :return: pandas dataframe
        """
        df = pd.read_csv(path_csv_file, header = 0)
        if self.downsampling == "Balanced":
            if self.cellnumber < df.shape[0]:
                df = df.sample(n = int(self.cellnumber), random_state = 42,
                               ignore_index = True)
            else:
                # self.log.info("It was not possible to downsample because the fixed number of events is greater than the original number of events. Decrease the threshold for downsampling.")
                pass
        barcode = []
        names = os.path.basename(path_csv_file)
        num_lines = df.shape[0]
        for _ in range(1, num_lines + 1):
            barcode.append("_".join([names.split(".")[0], str(_)]))
        df.index = barcode
        return df

    def concatenate_dataframe(self, info_file, csv_list):
        """

        :param csv_list:
        :return:
        """
        self.log.info("Files concatenation")
        # create empy list for save several df
        pandas_df_list = []
        # create list with anndata object
        # loop over csv file name
        for i in range(len(csv_list)):
            if self.filetype == "CSV":
                # append df in pandas_df_list
                pandas_df_list.append(self.create_df(csv_list[i]))
            else:
                pandas_df_list.append(self.read_fcs(csv_list[i]))
        # check header
        if all([len(pandas_df_list[0].columns.intersection(df.columns)) == pandas_df_list[0].shape[1]
                for df in pandas_df_list]):

            try:
                for i in range(len(pandas_df_list)):
                    # save column with Sample name in list
                    Sample_list = info_file["Sample"].tolist()
                    # check if Sample name are in the anndata index
                    if pandas_df_list[i].index[0][:-2] in Sample_list:
                        ann_tmp = anndata.AnnData(pandas_df_list[i])
                        ann_tmp.obs['Sample'] = pandas_df_list[i].index[0][:-2]
                        #
                        cell_type = info_file['Cell_type'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        ann_tmp.obs['Cell_type'] = ''.join(
                            e for e in cell_type.to_string().split(" ")[-1] if e.isalnum())
                        #
                        exp = info_file['EXP'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        ann_tmp.obs['EXP'] = ''.join(e for e in exp.to_string().split(" ")[-1] if e.isalnum())
                        #
                        id = info_file['ID'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        ann_tmp.obs['ID'] = ''.join(e for e in id.to_string().split(" ")[-1] if e.isalnum())
                        #
                        time_point = info_file['Time_point'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        # ann_tmp.obs['Time_point'] = time_point.to_string().split(" ")[-1]
                        ann_tmp.obs['Time_point'] = ''.join(
                            e for e in time_point.to_string().split(" ")[-1] if e.isalnum())
                        #

                        condition = info_file['Condition'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        ann_tmp.obs['Condition'] = ''.join(
                            e for e in condition.to_string().split(" ")[-1] if e.isalnum())
                        #
                        count = info_file['Count'].loc[info_file['Sample'] == pandas_df_list[i].index[0][:-2]]
                        ann_tmp.obs['Count'] = ''.join(e for e in count.to_string().split(" ")[-1] if e.isalnum())
                        self.anndata_list.append(ann_tmp)
                    else:
                        self.log.error(
                            "Error, this file {0} is not in the column Sample of Infofile.\nPlease check sample name and Infofile".format(
                                pandas_df_list[i].index[0][:-2]))
                        sys.exit(1)
                if len(self.anndata_list) == 1:
                    self.adata = self.anndata_list[0]
                    self.adata.layers['raw_value'] = self.adata.X
                else:
                    tmp = self.anndata_list[0]
                    self.anndata_list.pop(0)
                    self.adata = tmp.concatenate(self.anndata_list)
                    newheader = []
                    for _ in list(self.adata.var_names):
                        newheader.append(_.split(":: ")[-1])
                    self.adata.var_names = newheader
                    self.adata.layers['raw_value'] = self.adata.X
            except (ValueError, Exception):
                # print(traceback.format_exc())
                self.log.error("Error. Please check Info File Header or CSV header.")
                sys.exit(1)
        else:
            # print(traceback.format_exc())
            self.log.error("Error. Please check Info File Header or CSV header.")
            sys.exit(1)
        self.tmp_df = pd.DataFrame(self.adata.X, index = self.adata.obs.index)
        if self.downsampling == "Fixed":
            if self.cellnumber < self.adata.shape[0]:
                sc.pp.subsample(self.adata, n_obs = self.cellnumber, random_state = 42)
            else:
                pass
        self.cleaning.update({"Before QC":self.adata.shape[0]})
        self.log.info("{0} cells undergo to clustering analysis".format(self.adata.shape[0]))
        return self.adata

    def transformation(self):
        if (self.filetype == "FCS") and (self.arcsinh == True):
            self.adata.layers['arcin'] = scprep.transform.arcsinh(self.adata.X, cofactor=150)

    def create_barplot(self):
        """
        Create a barplot and export with the self.cleaning dictionary
        :return:
        """
        self.outfig = "/".join([self.output_folder, "".join(["Figures", self.tool])])
        self.createdir(self.outfig)
        self.QC_folder = "/".join([self.outfig, "QC_PLOTS"])
        self.createdir(self.QC_folder)
        ax = sns.barplot(data = pd.DataFrame.from_dict(self.cleaning, orient = 'index').reset_index(),
                    y = 0, x = 'index')
        ax.bar_label(ax.containers[0], fmt = '%.0f')
        plt.grid(False)
        plt.ylabel("Number of cells")
        plt.xlabel("Cleaning steps")
        plt.savefig(self.QC_folder+ "/cleaning.png", dpi = 300, bbox_inches = 'tight')

    def correct_scanorama(self):
        """
        This function runs Scanorama
        Returns: corrected adata
        """
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.adata_subset.layers['raw_value'] = self.adata_subset.X
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                         zero_center = True, copy = True).X
        self.anndata_list = [self.adata_subset[self.adata_subset.obs[self.batchcov] == i] for i in
                             self.adata_subset.obs[self.batchcov].unique()]
        self.corrected = scanorama.correct_scanpy(self.anndata_list,
                                                  return_dense = True,
                                                  return_dimred = True,
                                                  approx = False,
                                                  verbose = 0,
                                                  seed = 42)
        self.corrected_dataset = self.corrected[0].concatenate(self.corrected[1:],
                                                               join = 'inner',
                                                               batch_key = self.batchcov)
        self.corrected_dataset.layers['raw_value'] = self.adata_subset.X
        self.corrected_dataset.layers['scaled'] = self.adata_subset.layers['scaled']
        return self.corrected_dataset

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
            self.log.error("File does not exist.")
            sys.exit(1)

    def checkmarkers(self):
        """
        Check if marker in file is also a column of conc file
        :return:
        """
        # read marker file
        self.marker_array = [line.rstrip() for line in open(self.marker_list)]
        newmarker = []
        # if self.filetype == "CSV":
        for _ in self.marker_array:
            newmarker.append(_.split(":: ")[-1])
        self.marker_array = newmarker
        if len(self.anndata_list) > 1:
            # read concatenate file
            for i in range(len(self.marker_array)):
                if self.marker_array[i] in self.adata.var_names.to_list():
                    continue
                else:
                    self.log.error("Marker {} not found in Matrix.".format(self.marker_array[i]))
                    sys.exit(1)
            return self.marker_array
        else:
            return self.marker_array

    def splitmarker(self):
        """
        function for split marker in two list
        return: list of marker to include in analysis
        """
        self.marker = self.adata.var_names.to_list()
        self.markertoinclude = [i for i in self.marker if i not in self.marker_array]
        return self.markertoinclude

    def runumap(self):
        """
        Function for UMAP generation
        return: UMAP embedding
        """
        self.log.warning("PART 3")
        if self.runtime != 'Clustering':
            self.log.info("UMAP (Uniform Manifold Approximation and Projection) generation")
            reducer = umap.UMAP(random_state = 42, n_neighbors = self.n_neighbors, min_dist = self.mindist,
                                spread = self.spread)
            embedding = reducer.fit_transform(self.adata_subset.X)
            return embedding
        else:
            self.log.info("Skipping UMAP (Uniform Manifold Approximation and Projection) generation")

    def subsample_adata_plotting(self):
        """
        Function for subsample adata for plotting
        Returns: anndata object

        """
        if self.runtime != 'Clustering':
            adatas = [self.adata_subset[self.adata_subset.obs['pheno_leiden'].astype("category").isin([clust])] for clust in
                      self.adata_subset.obs['pheno_leiden'].astype("category").cat.categories]
            for dat in adatas:
                    sc.pp.subsample(dat, fraction = self.target_cells, random_state = 42)

            self.adata_downsampled = adatas[0].concatenate(*adatas[1:])
            # return self.adata_downsampled
        else:
            self.log.info("Skipping subsampling adata for plotting")

    def plot_umap(self):
        """
        Function per generation of pdf files with umap plot
        return: pdf files
        """
        if self.runtime == 'Full':
            # create output directory
            if self.tool != "FlowSOM":
                self.UMAP_folder = "/".join([self.outfig, "UMAP"])
                self.createdir(self.UMAP_folder)
                sc.settings.figdir = self.UMAP_folder
            else:
                sc.settings.figdir = self.UMAP_folder
            # set palette
            if len(self.adata_subset.obs["pheno_leiden"].unique()) < 28:
                self.palette = self.palette28
            else:
                self.palette = self.palette102
            # plot umap + clustering
            sc.pl.umap(self.adata_subset, color = "pheno_leiden",
                       legend_fontoutline = 2, show = False, add_outline = False, frameon = False,
                       title = "UMAP Plot", palette = self.palette,
                       s = 50, save = ".".join(["".join([str(self.tool), "_cluster"]), "pdf"]))
            sc.pl.umap(self.adata_subset, color = "pheno_leiden",
                       legend_fontoutline = 4, show = False, add_outline = False, frameon = False,
                       legend_loc = 'on data', title = "UMAP Plot", palette = self.palette,
                       s = 50, save = "_legend_on_data.".join(["".join([str(self.tool), "_cluster"]), self.fileformat]))
            # format svg
            sc.pl.umap(self.adata_subset, color = "pheno_leiden",
                       legend_fontoutline = 4, show = False, add_outline = False, frameon = False,
                       legend_loc = 'on data', title = "UMAP Plot", palette = self.palette,
                       s = 50, save = "_legend_on_data.".join(["".join([str(self.tool), "_cluster"]), 'svg']))
            # plot umap with info file condition
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    sc.pl.umap(self.adata_subset, color = _, legend_fontoutline = 2, show = False, add_outline = False,
                               frameon = False,
                               title = "UMAP Plot",
                               s = 50, save = ".".join(["_".join([str(self.tool), _]), "pdf"]))
                else:
                    continue
            # plot umap grouped with gray background
            for _ in ['Cell_type', 'EXP', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    for batch in list(self.adata_subset.obs[_].unique()):
                        sc.pl.umap(self.adata_subset, color = _, groups = [batch], na_in_legend = False,
                                   title = "UMAP Plot",
                                   legend_fontoutline = 2, show = False, add_outline = False, frameon = False,
                                   s = 50, save = ".".join(["_".join([_ + str(batch), _]), "pdf"]))
                else:
                    continue
            # scale data
            self.scaler = MinMaxScaler(feature_range = (0, 1))
            self.adata_subset.layers['scaled01'] = self.scaler.fit_transform(self.adata_subset.layers['raw_value'])
            for _ in list(self.adata_subset.var_names.unique()):
                if self.scaler is True:
                    sc.pl.umap(self.adata_subset, color = _, show = False, layer = "raw_value",
                               legend_fontoutline = 1, na_in_legend = False, s = 30,
                               title = _, cmap = 'turbo', groups = [_],
                               save = ".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                               )
                else:
                    sc.pl.umap(self.adata_subset, color = _, show = False, layer = "scaled01",
                               legend_fontoutline = 1, na_in_legend = False, s = 30,
                               title = _, cmap = 'turbo', groups = [_],
                               save = ".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                               )
        elif self.runtime == 'UMAP':
            self.outfig = "/".join([self.output_folder, "".join(["Figures", self.tool])])
            sc.settings.figdir = self.outfig
            scaler = MinMaxScaler(feature_range = (0, 1))
            self.adata_subset.layers['scaled01'] = scaler.fit_transform(self.adata_subset.layers['raw_value'])
            for _ in list(self.adata_subset.var_names.unique()):
                sc.pl.umap(self.adata_subset, color = _, show = False, layer = "scaled01",
                           legend_fontoutline = 1, na_in_legend = False, s = 30, frameon = False,
                           title = _, cmap = 'turbo', groups = [_],
                           save = ".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                           )
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    sc.pl.umap(self.adata_subset, color = _,
                               cmap = self.palette, legend_fontoutline = 2, show = False, add_outline = False,
                               frameon = False,
                               title = "UMAP Plot",
                               s = 50, save = ".".join(["_".join([str(self.tool), _]), "pdf"]))
                else:
                    continue
            for _ in ['Cell_type', 'EXP', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    for batch in list(self.adata_subset.obs[_].unique()):
                        sc.pl.umap(self.adata_subset, color = _, groups = [batch], na_in_legend = False,
                                   title = "UMAP Plot",
                                   legend_fontoutline = 2, show = False, add_outline = False, frameon = False,
                                   s = 50, save = ".".join(["_".join([_ + str(batch), _]), "pdf"])
                                   )
                else:
                    continue
            # self.plot_cell_obs()
        elif self.runtime == 'Clustering':
            pass

    def plot_cell_clusters(self):
        if self.runtime == 'Full':
            self.umap = pd.DataFrame(self.adata_downsampled.obsm['X_umap'], index = self.adata_downsampled.obs_names)
            clusters = self.adata_downsampled.obs['pheno_leiden']
            tsne = self.umap.copy()
            tsne.columns = ['x', 'y']

            # Cluster colors
            n_clusters = len(set(clusters))
            cluster_colors = pd.Series(
                sns.color_palette(self.palette, n_clusters), index = set(clusters)
            )

            # Set up figure
            n_cols = 6
            n_rows = int(np.ceil(n_clusters / n_cols))
            fig = plt.figure(figsize = [2 * n_cols, 2 * (n_rows + 2)], dpi = 300)
            gs = plt.GridSpec(
                n_rows + 2, n_cols, height_ratios = np.append([0.75, 0.75], np.repeat(1, n_rows))
            )

            # Clusters
            ax = plt.subplot(gs[0:2, 2:4])
            ax.scatter(tsne["x"], tsne["y"], s = 6, color = cluster_colors[clusters[tsne.index]])
            ax.set_axis_off()

            # Branch probabilities
            for i, cluster in enumerate(set(clusters)):
                row = int(np.floor(i / n_cols))
                ax = plt.subplot(gs[row + 2, i % n_cols])
                ax.scatter(tsne.loc[:, "x"], tsne.loc[:, "y"], s = 3, color = "lightgrey")
                cells = clusters.index[clusters == cluster]
                ax.scatter(
                    tsne.loc[cells, "x"],
                    tsne.loc[cells, "y"],
                    s = 3,
                    color = cluster_colors[cluster],
                )
                ax.set_axis_off()
                ax.set_title(cluster, fontsize = 10)
            fig.tight_layout()
            fig.savefig("".join([self.UMAP_folder, ".".join(["/umapCELL_clusters_all", self.fileformat])]))
            plt.close(fig)
        else:
            pass

    def plot_umap_expression(self):
        if self.runtime == 'Full':
            self.subsample_adata_plotting()
            self.adata_downsampled.obs['Clustering'] = self.adata_downsampled.obs['pheno_leiden'].astype(str)
            sc.pl.umap(self.adata_downsampled,
                       color = ['Clustering'] + list(self.adata_downsampled.var_names),
                       show = False,
                       layer = "scaled01",
                       legend_fontoutline = 1,frameon=False,
                       na_in_legend = False, s = 50, cmap = 'turbo',
                       save = ".".join(["".join([str(self.tool), "_ALL"]), self.fileformat])
                       )
            sc.pl.umap(self.adata_downsampled,
                       color = ['Clustering'] + list(self.adata_downsampled.var_names),
                       show = False,
                       layer = "scaled01",
                       legend_fontoutline = 1,frameon=False,
                       na_in_legend = False, s = 50, cmap = 'turbo',
                       save = ".".join(["".join([str(self.tool), "_ALL"]), 'svg'])
                       )
        else:
            pass
    def find_obs_not_unique(self):
        """
        Find obs columns that are not unique
        Returns: list of obs columns that are not unique

        """
        self.obs_not_unique = []
        for _ in self.adata_downsampled.obs.columns:
            if len(self.adata_downsampled.obs[_].unique()) == 1:
                self.obs_not_unique.append(_)
            else:
                continue
        return self.obs_not_unique

    def plot_cell_obs(self):
        if self.runtime != 'Clustering':
            for _ in ['Cell_type', 'EXP', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    sc.pl.umap(self.adata_downsampled,
                               color=['Clustering',_],
                               show=False,
                               layer="scaled01",
                               legend_fontoutline=1, frameon=False,
                               na_in_legend=False, s=50, cmap='turbo',
                               save=".".join(["".join([str(self.tool), _+"_ALL"]), self.fileformat])
                               )
                else:
                    continue
        else:
            pass

    def matrixplot(self):
        """
        Function for the generation of matrixplot sc.pl.matrixplot
        return: matrixplot
        """
        self.matrixplot_folder = "/".join([self.outfig, "HEATMAP"])
        self.createdir(self.matrixplot_folder)
        sc.settings.figdir = self.matrixplot_folder
        if self.runtime != 'UMAP':
            sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                             dendrogram = True, vmin = -2, vmax = 2, cmap = 'RdBu_r', layer = "scaled",
                             show = False, swap_axes = False, return_fig = False,
                             save = ".".join(["matrixplot_mean_z_score", self.fileformat]))
            sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                             dendrogram = True, vmin = -2, vmax = 2, cmap = 'RdBu_r', layer = "scaled",
                             show = False, swap_axes = False, return_fig = False,
                             save = ".".join(["matrixplot_mean_z_score", 'svg']))
            sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                             dendrogram = True, cmap = 'Blues', standard_scale = 'var',
                             colorbar_title = 'column scaled\nexpression', layer = "scaled",
                             swap_axes = False, return_fig = False,
                             show = False,
                             save = ".".join(["matrixplot_column_scaled_expression", self.fileformat]))
        else:
            pass

    def createfcs(self):
        """
        Function for the generation of fcs files and FlowAI QC
        return: Anndata with filtered cells
        """
        try:
            self.log.info("Perform Flow Auto QC with FlowAI tool.")
            df = pd.DataFrame(self.adata.X, columns = self.adata.var.index, index = self.adata.obs.index)
            self.concatenate_fcs = "/".join([self.output_folder,
                                             "Test_ConcatenatedCells.fcs"])
            if 'Time' in df.columns:
                fcsy.write_fcs(path = self.concatenate_fcs, df = df)
                subprocess.check_call(['Rscript', '--vanilla',
                                       self.path_flowai, self.concatenate_fcs,
                                       self.output_folder], stdout = self.fnull, stderr = self.fnull)
                df = fcsy.read_fcs("/".join([self.output_folder,
                                            "Test_ConcatenatedCells_concatenate_after_QC.fcs"]))
                df.set_index(self.adata.obs.index, inplace = True)
                self.adata.obs['remove_from_FM'] = df['remove_from_FM']
                self.adata = self.adata[(self.adata.obs['remove_from_FM'] < 10000), :]
                self.adata.layers['raw_value'] = self.adata.X
                self.log.info("{0} cells after FlowAI analysis".format(self.adata.shape[0]))
                self.cleaning.update({ "After QC": self.adata.shape[0]})
                self.create_barplot()
                return self.adata
            else:
                # self.log.info("Time channel not found. Skip QC")
                pass
        except:
            # print(traceback.format_exc())
            # self.log.info("Impossible to complete Flow Auto QC. Check Time channel.")
            pass

    def runphenograph(self):
        """
        Function for execution of phenograph analysis
        :return:
        """
        self.log.warning("PART 2")
        self.log.info("Phenograph Clustering")
        self.createfcs()
        self.log.info("Markers used for Phenograph clustering:")
        for i in self.markertoinclude:
            self.log.info(" + " + i)
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        if len(self.marker_array):
            self.log.info("Markers excluded for Phenograph clustering:")
        for i in self.marker_array:
            self.log.info(" - " + i)
        if (self.filetype == "FCS") and (self.arcsinh == True):
            self.adata_subset = scprep.transform.arcsinh(self.adata_subset.X, cofactor=150)
        else:
            pass
        if self.runtime != 'Clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        self.adata_subset.obs['pheno_leiden'], graph, Q = pg.cluster(self.adata_subset.X, k = int(self.k_coef),
                                                                     seed = 42,
                                                                     clustering_algo = "leiden",
                                                                     directed = True, primary_metric = "euclidean",
                                                                     q_tol = 0.05,
                                                                     prune = False, min_cluster_size = 1,
                                                                     n_jobs = int(self.thread))
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden'].values
        self.adata.obs['Phenograph_cluster'] = self.adata_subset.obs['pheno_leiden'].values
        if self.runtime != 'Clustering':
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        self.plot_umap_expression()
        self.plot_frequency()
        self.plot_cell_clusters()
        # self.plot_cell_obs()
        self.matrixplot()
        return self.adata

    def runvia(self):
        """
        function for execution of
        :return:
        """
        self.log.warning("PART 2")
        self.log.info("VIA Clustering")
        self.createfcs()
        self.log.info("Markers used for VIA clustering:")
        for i in self.markertoinclude:
            self.log.info(" + " + i)
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        if len(self.marker_array):
            self.log.info("Markers excluded for VIA clustering:")
        for i in self.marker_array:
            self.log.info(" - " + i)
        if (self.filetype == "FCS") and (self.arcsinh == True):
            self.adata_subset = scprep.transform.arcsinh(self.adata_subset.X, cofactor=150)
        else:
            pass
        if self.runtime != 'Clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        p = via.VIA(self.adata_subset.X, random_seed = 42, knn = int(self.knn), root_user = self.root_user,
                    jac_std_global = 'median', jac_weighted_edges = False, distance = 'l2',
                    small_pop = 10, resolution_parameter = self.resolution,
                    num_threads = int(self.thread))
        p.run_VIA()
        self.adata_subset.obs['pheno_leiden'] = [str(i) for i in p.labels]
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden'].values
        self.adata.obs['VIA_cluster'] = self.adata_subset.obs['pheno_leiden'].values
        if self.runtime != 'Clustering':
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        self.plot_umap_expression()
        self.plot_frequency()
        self.plot_cell_clusters()
        # self.plot_cell_obs()
        self.matrixplot()
        return self.adata

    def runflowsom(self):
        """
        function for execution of
        :return:
        """
        self.log.warning("PART 2")
        self.log.info("FlowSOM Clustering")
        self.createfcs()
        self.log.info("Markers used for FlowSOM clustering:")
        for i in self.markertoinclude:
            self.log.info(" + " + i)
        if len(self.marker_array):
            self.log.info("Markers excluded for FlowSOM clustering:")
        for i in self.marker_array:
            self.log.info(" - " + i)
        if (self.filetype == "FCS") and (self.arcsinh == True):
            self.adata_subset = scprep.transform.arcsinh(self.adata_subset.X, cofactor=150)
        else:
            pass
        if self.runtime != 'Clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset = self.adata[:, self.markertoinclude].copy()
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset = self.adata[:, self.markertoinclude].copy()
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                                 zero_center = True, copy = True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        self.adata_subset.X = self.adata_subset.layers['raw_value']
        ###
        self.tmp = self.adata_subset.to_df()
        self.tmp = self.tmp.astype(int)
        self.tmp.to_csv(self.output_folder+"/tmp.csv", header=True,
                        index=True, sep=',', mode='w')
        self.UMAP_folder = "/".join([self.outfig, "UMAP"])
        self.createdir(self.UMAP_folder)
        subprocess.check_call(['Rscript', '--vanilla',
                               self.path_flowsom, self.output_folder+"/tmp.csv",
                               self.output_folder,self.UMAP_folder,self.maxclus], stdout=self.fnull, stderr=self.fnull)
        self.flowsomDF = pd.read_csv(self.output_folder+"/output_flowsom.csv", sep=',', header=0, index_col=0)
        self.adata_subset.obs['Clusters'] = self.flowsomDF['Clusters'].values
        self.adata_subset.obs['Metaclusters'] = self.flowsomDF['Metaclusters'].values
        self.adata.obs['Cluster_Flowsom'] = self.adata_subset.obs['Clusters'].astype('category')
        self.adata.obs['MetaCluster_Flowsom'] = self.adata_subset.obs['Metaclusters'].astype('category')
        self.adata_subset.obs['pheno_leiden'] = self.flowsomDF['Metaclusters'].values
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype("category")
        self.adata.obs['cluster'] =self.flowsomDF['Metaclusters'].values
        self.adata_subset.X = self.adata_subset.layers['scaled']
        if self.runtime == 'Full':
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        self.plot_umap_expression()
        self.plot_frequency()
        self.plot_cell_clusters()
        # self.plot_cell_obs()
        self.matrixplot()
        return self.adata

    def generation_concatenate(self):
        """
        Function to concatenate the results of the clustering and the original adata object
        return: adata object with the results of the clustering
        """
        if self.runtime == 'Full':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns = self.adata.var_names,
                                                index = self.adata.obs.index).astype(int),
                                   pd.DataFrame(self.adata.obsm['X_umap'], columns = ['UMAP_1', 'UMAP_2'],
                                                index = self.adata.obs.index),
                                   right_index = True,
                                   left_index = True)
            pd.merge(self.tmp_df, self.adata.obs[['cluster',
                                                  'Sample', 'Cell_type',
                                                  'EXP',
                                                  'ID', 'Time_point',
                                                  'Condition']], left_index = True,
                     right_index = True).to_csv("/".join([self.output_folder, ".".join([self.analysis_name, "csv"])]),
                                                header = True, index = False)
        elif self.runtime == 'UMAP':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns = self.adata.var_names,
                                                index = self.adata.obs.index).astype(int),
                                   pd.DataFrame(self.embedding, columns = ['UMAP_1', 'UMAP_2'],
                                                index = self.adata.obs.index),
                                   right_index = True,
                                   left_index = True)
            pd.merge(self.tmp_df, self.adata.obs[['Sample', 'Cell_type',
                                                  'EXP',
                                                  'ID', 'Time_point',
                                                  'Condition']], left_index = True,
                     right_index = True).to_csv("/".join([self.output_folder, ".".join([self.analysis_name, "csv"])]),
                                                header = True, index = False)
        elif self.runtime == 'Clustering':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns = self.adata.var_names,
                                                index = self.adata.obs.index).astype(int),
                                   self.adata.obs[['cluster', 'Sample', 'Cell_type',
                                                   'EXP',
                                                   'ID', 'Time_point',
                                                   'Condition']],
                                   right_index = True,
                                   left_index = True)
            self.tmp_df.to_csv("/".join([self.output_folder, ".".join([self.analysis_name, "csv"])]), header = True,
                               index = False)

    def plot_frequency(self):
        """
        Plot barplot with frequency
        return: barplot with frequency
        """
        if self.runtime != 'UMAP':
            self.FREQUENCY_folder = "/".join([self.outfig, "BARPLOT_FREQUENCY"])
            self.createdir(self.FREQUENCY_folder)
            fig, (ax1) = plt.subplots(1, 1, figsize = (17 / 2.54, 17 / 2.54))
            ax1 = self.adata_subset.obs.groupby("pheno_leiden")["Sample"].value_counts(
                normalize = True).unstack().plot.barh(
                stacked = True,
                legend = False,
                ax = ax1,
                color = self.palette)
            ax1.set_xlabel("Percentage Frequency")
            ax1.set_ylabel("Cluster")
            ax1.grid(False)
            ax1.legend(bbox_to_anchor = (1.2, 1.0))
            ### save
            fig.savefig("/".join([self.FREQUENCY_folder , ".".join(["ClusterFrequencyNormalized", self.fileformat])]),
                        dpi = self.dpi, bbox_inches = 'tight',
                        format = self.fileformat)
            fig.savefig("/".join([self.FREQUENCY_folder , ".".join(["ClusterFrequencyNormalized", 'svg'])]),
                        dpi = self.dpi, bbox_inches = 'tight',
                        format = 'svg')
            #
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    fig, (ax3) = plt.subplots(1, 1, figsize = (17 / 2.54, 17 / 2.54))
                    ax3 = self.adata_subset.T.var.groupby(_)["pheno_leiden"].value_counts(
                        normalize = True).unstack().plot.barh(
                        stacked = True,
                        legend = False,
                        color = self.palette,
                        ax = ax3,
                        fontsize = 5)
                    ax3.set_xlabel("Cluster Percentage Frequency")
                    ax3.set_ylabel(_)
                    ax3.grid(False)
                    ax3.legend(bbox_to_anchor = (1.2, 1.0))
                    fig.savefig("/".join(
                        [self.FREQUENCY_folder , ".".join(["".join([_, "ClusterFrequencyNormalized"]), self.fileformat])]),
                                dpi = self.dpi, bbox_inches = 'tight',
                                format = self.fileformat)
            #
            fig, (ax2) = plt.subplots(1, 1, figsize = (17 / 2.54, 17 / 2.54))
            ax2 = self.adata_subset.obs.groupby("pheno_leiden")["Sample"].value_counts(
                normalize = False).unstack().plot.barh(
                stacked = True,
                legend = False,
                ax = ax2,
                color = self.palette)
            ax2.set_xlabel("Relative Frequency")
            ax2.set_ylabel("Cluster")
            ax2.grid(False)
            ax2.legend(bbox_to_anchor = (1.2, 1.0))
            fig.savefig("/".join([self.FREQUENCY_folder , ".".join(["ClusterFrequencyNotNormalized", self.fileformat])]),
                        dpi = self.dpi, bbox_inches = 'tight',
                        format = self.fileformat)
            fig.savefig("/".join([self.FREQUENCY_folder , ".".join(["ClusterFrequencyNotNormalized", 'svg'])]),
                        dpi = self.dpi, bbox_inches = 'tight',
                        format = 'svg')
            self.plot_frequency_ptz()
        else:
            pass

    def plot_frequency_ptz(self):
        """
        Function to plot frequency of clusters per time point and per condition
        Returns: barplot with frequency per time point and per condition
        """
        if len(self.adata_subset.obs["ID"].unique()) > 1:
            self.dfxlinkage = self.adata_subset.obs.groupby("ID")["pheno_leiden"].value_counts(
                normalize = True).unstack() * 100
            self.dfxlinkage.fillna(0, inplace = True)
            Z = linkage(self.dfxlinkage, 'ward',
                        optimal_ordering = True)
            dn = dendrogram(Z, get_leaves = True, orientation = 'left', labels = self.dfxlinkage.index,
                            no_plot = True)
            fig, (ax1, ax2) = plt.subplots(1, 2, constrained_layout = True, figsize = (20, 10))
            dn = dendrogram(Z, get_leaves = True, orientation = 'left', labels = self.dfxlinkage.index,
                            color_threshold = 0, above_threshold_color = 'k', ax = ax1)
            self.dfxlinkage.loc[dn['ivl']].plot.barh(legend = False, stacked = True, ax = ax2, color = self.palette)
            ax1.set(yticklabels = [])
            ax1.set(xticklabels = [])
            ax1.grid(False)
            ax2.tick_params(left = False)
            ax2.grid(False)
            ax1.axis('off')
            ax2.set_ylabel(" ")
            ax2.set_xlabel("Percentage Frequency")
            ax2.legend(bbox_to_anchor = (1.2, 1.0), title = 'Cluster')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            fig.savefig("/".join([self.FREQUENCY_folder , ".".join(["SampleFrequencyClusterized", self.fileformat])]),
                        dpi = self.dpi, bbox_inches = 'tight',
                        format = self.fileformat)
        else:
            pass

    def createdir(self, dirpath):
        """
        Make dir function and check if directory is already exists
        :param dirpath: string with path and directory name
        :return:
        """
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)

    def groupbycluster(self):
        """
        Function for generation of csv with different clusters
        :return:
        """
        # make dir
        self.createdir("/".join([self.output_folder, "".join(["CSVcluster", self.tool])]))
        self.createdir("/".join([self.output_folder, "".join(["FCScluster", self.tool])]))
        self.adata.obs['cluster'] = self.adata.obs['cluster'].astype(int)
        for _ in range(self.adata.obs['cluster'].unique().min(), self.adata.obs['cluster'].unique().max() + 1):
            self.tmp = self.adata[self.adata.obs['cluster'].isin([_])].to_df()
            self.tmp = self.tmp.astype(int)
            if self.runtime == 'Full':
                self.tmp['UMAP_1'] = pd.DataFrame(
                    self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:, 0]).set_index(self.tmp.index)
                self.tmp['UMAP_2'] = pd.DataFrame(
                    self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:, 1]).set_index(self.tmp.index)
            if (self.tool == "Phenograph"):
                self.tmp['Phenograph'] = _
            elif (self.tool == "VIA"):
                self.tmp['VIA'] = _
            else:
                self.tmp['FlowSOM'] = _
                self.tmp['MetaCluster_FlowSOM'] = self.adata[self.adata.obs['cluster'].isin([_])].obs[
                    'Cluster_Flowsom'].values
            self.tmp.to_csv("/".join([self.output_folder, "".join(["CSVcluster", self.tool]),
                                      "".join([self.analysis_name, "_", str(_), ".csv"])]), header = True,
                            index = False)
            fcsy.write_fcs(self.tmp, "/".join([self.output_folder, "".join(["FCScluster", self.tool]),
                                               "".join([self.analysis_name, "_", str(_), ".fcs"])]))

    def groupbysample(self):
        """
        Function for generation of csv with different clusters
        """
        # make dir
        self.createdir("/".join([self.output_folder, "".join(["CSVsample", self.tool])]))
        self.createdir("/".join([self.output_folder, "".join(["FCSsample", self.tool])]))
        # tmp df
        self.tmp = self.adata.to_df()
        self.tmp = self.tmp.astype(int)
        # change index
        self.tmp.set_index(self.adata.obs['Sample'], inplace = True)
        if self.runtime != 'UMAP':
            self.adata.obs['cluster'] = self.adata.obs['cluster'].astype(int)
            self.createdir("/".join([self.output_folder, "".join(["ClusterFrequency", self.tool])]))
        # create columns
        if self.runtime != 'Clustering':
            self.tmp['UMAP_1'] = self.adata.obsm['X_umap'][:, 0]
            self.tmp['UMAP_2'] = self.adata.obsm['X_umap'][:, 1]
        if self.runtime != 'UMAP':
            self.tmp[self.tool] = self.adata.obs['cluster'].values
            if self.tool == 'FlowSOM':
                self.tmp['MetaCluster_FlowSOM'] = self.adata.obs['Cluster_Flowsom'].values
            else:
                pass
            self.tmp["cluster"] = self.adata.obs['cluster'].values
            # get unique filenames
            unique_filename = self.adata.obs['Sample'].unique()
            # get unique number of cluster
            unique_Phenograph = self.adata.obs['cluster'].unique()
            #
            dfCounts = pd.DataFrame(index = range(min(unique_Phenograph), max(unique_Phenograph) + 1))
            # generate Tot_percentage file
            for i in range(len(unique_filename)):
                dfCounts[unique_filename[i]] = self.tmp.loc[self.tmp.index == unique_filename[i]].cluster.value_counts(
                    normalize = True).reindex(self.tmp.cluster.unique(), fill_value = 0)
            # compute percentage
            dfCounts = dfCounts * 100
            # save
            dfCounts.index.name = 'Cluster'
            dfCounts.to_csv("/".join(["/".join([self.output_folder, "".join(["ClusterFrequency", self.tool])]),
                                      "".join(["TotalPercentage", ".csv"])]))
            # create empty dataframe
            dfCounts = pd.DataFrame(index = range(min(unique_Phenograph), max(unique_Phenograph) + 1))
            # generate Tot_counts file
            for i in range(len(unique_filename)):
                dfCounts[unique_filename[i]] = self.tmp.loc[
                    self.tmp.index == unique_filename[i]].cluster.value_counts().reindex(
                    self.tmp.cluster.unique(), fill_value = 0)
            # save
            dfCounts.index.name = 'Cluster'
            dfCounts.to_csv("/".join(["/".join([self.output_folder, "".join(["ClusterFrequency", self.tool])]),
                                      "".join(["TotalCounts", ".csv"])]))
            del self.tmp['cluster']
            # save samples
            for i in range(len(unique_filename)):
                dfCounts[unique_filename[i]] = self.tmp.loc[self.tmp.index == unique_filename[i]].to_csv(
                    "/".join([self.output_folder, "".join(["CSVsample", self.tool]),
                              "".join([str(unique_filename[i]), "_", self.analysis_name,
                                       ".csv"])]),
                    header = True, index = False)
                fcsy.write_fcs(self.tmp.loc[self.tmp.index == unique_filename[i]],
                               "/".join([self.output_folder, "".join(["FCSsample", self.tool]),
                                         "".join([str(unique_filename[i]), "_", self.analysis_name,
                                                  ".fcs"])])
                               )
        else:
            unique_filename = self.adata.obs['Sample'].unique()
            for i in range(len(unique_filename)):
                self.tmp.loc[self.tmp.index == unique_filename[i]].to_csv(
                    "/".join([self.output_folder, "".join(["CSVsample", self.tool]),
                              "".join([str(unique_filename[i]), "_", self.analysis_name,
                                       ".csv"])]),
                    header = True, index = False)
                fcsy.write_fcs(self.tmp.loc[self.tmp.index == unique_filename[i]],
                               "/".join([self.output_folder, "".join(["FCSsample", self.tool]),
                                         "".join([str(unique_filename[i]), "_", self.analysis_name, ".fcs"])]))

    def runtimeumap(self):
        """
        Function for execution of phenograph analysis
        :return:
        """
        self.log.warning("PART 2")
        self.log.info("UMAP Dimensional Reduction")

        self.log.info("Markers used for UMAP computation:")
        for i in self.markertoinclude:
            self.log.info(" + " + i)
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        if len(self.marker_array):
            self.log.info("Markers excluded for UMAP computation:")
        for i in self.marker_array:
            self.log.info(" - " + i)
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value = 6,
                                                         zero_center = True, copy = True).X
        # self.adata_subset.X = self.adata_subset.layers['scaled']
        if self.scanorama is True:
            self.adata_subset = self.correct_scanorama()
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        else:
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        return self.adata

    def exporting(self):
        """
        Export to h5ad file.
        """
        self.log.warning("PART 4")
        self.log.info("Output Generation")
        self.scaler = MinMaxScaler(feature_range = (0, 1))
        try:
            del self.adata.obs['remove_from_FM']
        except:
            pass
        if self.runtime != 'UMAP':
            if self.tool == "Phenograph":
                self.adata.obs[self.tool + "_" + str(self.k_coef)] = self.adata.obs['cluster'].astype("str")
                del self.adata.obs['cluster']
                del self.adata.obs[self.tool + "_" + str(self.k_coef)]
                self.adata.obs["".join([str(self.tool), "_cluster"])] = self.adata.obs[
                    "".join([str(self.tool), "_cluster"])].astype('category')
                self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
                self.adata.X = self.adata.layers['scaled01']
                self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
                self.adata.write("/".join([self.output_folder, ".".join([self.analysis_name, "h5ad"])]))
            elif self.tool == "VIA":
                self.adata.obs[self.tool + "_" + str(self.knn)] = self.adata.obs['cluster'].astype("str")
                del self.adata.obs['cluster']
                del self.adata.obs[self.tool + "_" + str(self.knn)]
                self.adata.obs["".join([str(self.tool), "_cluster"])] = self.adata.obs[
                    "".join([str(self.tool), "_cluster"])].astype('category')
                self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
                self.adata.X = self.adata.layers['scaled01']
                self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
                self.adata.write("/".join([self.output_folder, ".".join([self.analysis_name, "h5ad"])]))
            else:
                del self.adata.obs['cluster']
                self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
                self.adata.X = self.adata.layers['scaled01']
                self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
                self.adata.write("/".join([self.output_folder, ".".join([self.analysis_name, "h5ad"])]))
                try:
                    os.remove(self.output_folder+ "tmp.csv")
                    os.remove(self.output_folder + "output_flowsom.csv")
                except:
                    pass
        else:
            self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
            self.adata.X = self.adata.layers['scaled01']
            self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
            self.adata.write("/".join([self.output_folder, ".".join([self.analysis_name, "h5ad"])]))
        try:
            os.remove(self.output_folder + "".join([self.analysis_name, "_ConcatenatedCells_concatenate_after_QC.fcs"]))
            os.remove(self.output_folder + "".join([self.analysis_name, "_ConcatenatedCells.fcs"]))
        except:
            pass
        self.log.warning("PART 5")
