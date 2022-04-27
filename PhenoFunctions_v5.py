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
import matplotlib
import scanorama
import scipy
import matplotlib.pyplot as plt
import seaborn as sb
from sklearn.preprocessing import MinMaxScaler
matplotlib.use('Agg')


tmp = tempfile.NamedTemporaryFile()
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=300, facecolor='white', dpi_save=330,
                              figsize=(10, 10))
sc.settings.verbosity = 0
warnings.filterwarnings("ignore", category=FutureWarning)


class Cytophenograph:
    def __init__(self, info_file, input_folder, output_folder, k_coef, marker_list, analysis_name, thread, tool, batch,
                 batchcov, mindist, spread, runtime):
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
        self.listmarkerplot = None
        if 0.001 <= float(mindist) <= 1.0:
            self.mindist = float(mindist)
        else:
            self.log.error("Error. Min_dist parameter should be between 0.1 and 1.")
            sys.exit(1)
        if 0.001 <= float(spread) <= 1.0:
            self.spread = float(spread)
        else:
            self.log.error("Error. Min_dist parameter should be between 0.1 and 1.")
            sys.exit(1)
        format = logging.Formatter("%(asctime)s %(threadName)-11s %(levelname)-10s %(message)s")
        #
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        self.log.addHandler(ch)
        #
        fh = logging.FileHandler("/".join([self.output_folder, "log.txt"]), "w")
        fh.setFormatter(format)
        self.log.addHandler(fh)
        self.log.info("Name of this analysis: {}".format(marker_list))
        self.log.info("Input folder: {}".format(input_folder))
        self.log.info("Output folder: {}".format(output_folder))
        self.log.info("Info file: {}".format(info_file))
        self.log.info("Phenograph K-coef: {}".format(k_coef))
        self.log.info("Marker list file: {}".format(marker_list))
        self.log.info("Clustering tool option: {}".format(tool))
        self.log.info("Runtime option: {}".format(runtime))
        if self.scanorama is True:
            self.log.info("Covariate selected for batch correction is: {}".format(batchcov))
        self.log.info("UMAP min_dist is: {}".format(mindist))
        self.log.info("UMAP spread is: {}".format(spread))
        self.palette28 = ["#08519c",  # Blue 1
         "#ff7f0e",  # Orange 1
         "#1f6836",  # Green 1
         "#514888",  # Purple 1
         "#b30000",  # Red 1
         "#5a3730",  # Brown 1
         "#d638a6",  # Pink 1
         "#595959",  # Grey 1
         "#7c7c16",  # Green 1
         "#77e398",  # Light Green
         "#3182bd",  # Blue 2
         "#9e3a09",  # Orange 2
         "#31a354",  # Green 2
         "#756bb1",  # Purple 2
         "#ff0000",  # Red 2
         "#8c564b",  # Brown 2
         "#e377c2",  # Pink 2
         "#808080",  # Grey 2
         "#bcbd22",  # Green 2
         "#85b5d3",  # Blue 3
         "#ffa85b",  # Orange 3
         "#55cc79",  # Green 3
         "#a49dcb",  # Purple 3
         "#ff4d4d",  # Red 3
         "#b37c71",  # Brown 3
         "#f0b6de",  # Pink 3
         "#a6a6a6",  # Grey 3
         "#dedf4d"]  # Green 3
        self.palette102 = [
            "#FFFF00",
            "#1CE6FF",
            "#FF34FF",
            "#FF4A46",
            "#008941",
            "#006FA6",
            "#A30059",
            "#FFDBE5",
            "#7A4900",
            "#0000A6",
            "#63FFAC",
            "#B79762",
            "#004D43",
            "#8FB0FF",
            "#997D87",
            "#5A0007",
            "#809693",
            "#6A3A4C",
            "#1B4400",
            "#4FC601",
            "#3B5DFF",
            "#4A3B53",
            "#FF2F80",
            "#61615A",
            "#BA0900",
            "#6B7900",
            "#00C2A0",
            "#FFAA92",
            "#FF90C9",
            "#B903AA",
            "#D16100",
            "#DDEFFF",
            "#000035",
            "#7B4F4B",
            "#A1C299",
            "#300018",
            "#0AA6D8",
            "#013349",
            "#00846F",
            "#372101",
            "#FFB500",
            "#C2FFED",
            "#A079BF",
            "#CC0744",
            "#C0B9B2",
            "#C2FF99",
            "#001E09",
            "#00489C",
            "#6F0062",
            "#0CBD66",
            "#EEC3FF",
            "#456D75",
            "#B77B68",
            "#7A87A1",
            "#788D66",
            "#885578",
            "#FAD09F",
            "#FF8A9A",
            "#D157A0",
            "#BEC459",
            "#456648",
            "#0086ED",
            "#886F4C",
            "#34362D",
            "#B4A8BD",
            "#00A6AA",
            "#452C2C",
            "#636375",
            "#A3C8C9",
            "#FF913F",
            "#938A81",
            "#575329",
            "#00FECF",
            "#B05B6F",
            "#8CD0FF",
            "#3B9700",
            "#04F757",
            "#C8A1A1",
            "#1E6E00",
            "#7900D7",
            "#A77500",
            "#6367A9",
            "#A05837",
            "#6B002C",
            "#772600",
            "#D790FF",
            "#9B9700",
            "#549E79",
            "#FFF69F",
            "#201625",
            "#72418F",
            "#BC23FF",
            "#99ADC0",
            "#3A2465",
            "#922329",
            "#5B4534",
            "#FDE8DC",
            "#404E55",
            "#0089A3",
            "#CB7E98",
            "#A4E804",
            "#324E72",
        ]

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

    def concatenate_dataframe(self, info_file, csv_list):
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
                            "Error, this file {0} is not in the column Sample of Infofile. \n Please check sample name and Infofile".format(
                                pandas_df_list[i].index[0][:-2]))
                        sys.exit(1)
                if len(self.anndata_list) == 0:
                    self.adata = self.anndata_list[0]
                    self.adata.layers['raw_value'] = self.adata.X
                else:
                    tmp = self.anndata_list[0]
                    self.anndata_list.pop(0)
                    self.adata = tmp.concatenate(self.anndata_list)
                    self.adata.layers['raw_value'] = self.adata.X
            except (ValueError, Exception):
                self.log.error("Error. Please check Info File Header or CSV header.")
                sys.exit(1)
        else:
            self.log.error("Error. Please check Info File Header or CSV header.")
            sys.exit(1)
        self.tmp_df = pd.DataFrame(self.adata.X, index=self.adata.obs.index)
        if sum('::' in s for s in list(self.adata.var_names)) == len(list(self.adata.var_names)):
            for _ in range(len(list(self.adata.var_names))):
                self.newheader.append(list(self.adata.var_names)[_].split("::")[-1].upper())
            self.tmp_df.columns = [x.strip(' ') for x in self.newheader]
        else:
            self.tmp_df.columns = self.adata.var_names
        return self.adata

    def correct_scanorama(self):
        """
        This function runs Scanorama
        Returns: corrected adata
        """
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.adata_subset.layers['raw_value'] = self.adata_subset.X
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                         zero_center=True, copy=True).X
        self.anndata_list = [self.adata_subset[self.adata_subset.obs[self.batchcov] == i] for i in
                             self.adata_subset.obs[self.batchcov].unique()]
        self.corrected = scanorama.correct_scanpy(self.anndata_list,
                                                  return_dense=True,
                                                  return_dimred=True,
                                                  approx=False,
                                                  verbose=0,
                                                  seed=42)
        self.corrected_dataset = self.corrected[0].concatenate(self.corrected[1:],
                                                               join='inner',
                                                               batch_key=self.batchcov)
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
        if self.runtime != 'clustering':
            self.log.info("Part3: UMAP (Uniform Manifold Approximation and Projection) generation")
            reducer = umap.UMAP(random_state=42, n_neighbors=self.n_neighbors, min_dist=self.mindist, spread=1.0)
            embedding = reducer.fit_transform(self.adata_subset.X)
            return embedding
        else:
            self.log.info("Part3: skipping UMAP (Uniform Manifold Approximation and Projection) generation")

    def plot_umap(self):
        """
        Function per generation of pdf files with umap plot
        """
        if self.runtime == 'full':
            self.createdir("/".join([self.output_folder, "".join(["Figures", self.tool])]))
            self.outfig = "/".join([self.output_folder, "".join(["Figures", self.tool])])
            sc.settings.figdir = self.outfig
            # set palette
            if len(self.adata_subset.obs["pheno_leiden"].unique()) < 28:
                self.palette = self.palette28
            else:
                self.palette = self.palette102
            # umap cluster
            sc.pl.umap(self.adata_subset, color="pheno_leiden",
                       legend_fontoutline=2, show=False, add_outline=False, frameon=False,
                       title="UMAP Plot",palette=self.palette,
                       s=50, save=".".join(["".join([str(self.tool), "_cluster"]), "pdf"]))
            sc.pl.umap(self.adata_subset, color="pheno_leiden",
                       legend_fontoutline=2, show=False, add_outline=False, frameon=False,
                       legend_loc='on data', title="UMAP Plot",palette=self.palette,
                       s=50, save="_legend_on_data.".join(["".join([str(self.tool), "_cluster"]), self.fileformat]))
            # umap obs
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    sc.pl.umap(self.adata_subset, color=_, legend_fontoutline=2, show=False, add_outline=False,
                               frameon=False,
                               title="UMAP Plot",
                               s=50, save=".".join(["_".join([str(self.tool), _]), "pdf"]))
            # scale data
            self.scaler = MinMaxScaler(feature_range=(0, 1))
            self.adata_subset.layers['scaled01'] = self.scaler.fit_transform(self.adata_subset.layers['raw_value'])
            for _ in list(self.adata_subset.var_names.unique()):
                if self.scaler is True:
                    sc.pl.umap(self.adata_subset, color=_, show=False, layer="raw_value",
                               legend_fontoutline=1, na_in_legend=False, s=30,
                               title=_, cmap='turbo', groups=[_],
                               save=".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                               )
                else:
                    sc.pl.umap(self.adata_subset, color=_, show=False, layer="scaled01",
                               legend_fontoutline=1, na_in_legend=False, s=30,
                               title=_, cmap='turbo', groups=[_],
                               save=".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                               )
            self.listmarkerplot = ['pheno_leiden']
            for _ in range(len(list(self.adata_subset.var_names))):
                self.listmarkerplot.append(list(self.adata_subset.var_names)[_])
            sc.pl.umap(self.adata_subset, color=self.listmarkerplot, show=False, layer="scaled01",
                       legend_fontoutline=1, na_in_legend=False, s=30, cmap='turbo',
                       save=".".join(["".join([str(self.tool), "_all"]), self.fileformat])
                       )
        elif self.runtime == 'umap':
            self.createdir("/".join([self.output_folder, "".join(["Figures", "_runtimeUMAP"])]))
            self.outfig = "/".join([self.output_folder, "".join(["Figures", "_runtimeUMAP"])])
            sc.settings.figdir = self.outfig
            scaler = MinMaxScaler(feature_range=(0, 1))
            self.adata_subset.layers['scaled01'] = scaler.fit_transform(self.adata_subset.layers['raw_value'])
            for _ in list(self.adata_subset.var_names.unique()):
                sc.pl.umap(self.adata_subset, color=_, show=False, layer="scaled01",
                           legend_fontoutline=1, na_in_legend=False, s=30, frameon=False,
                           title=_, cmap='turbo', groups=[_],
                           save=".".join([''.join(e for e in _ if e.isalnum()), self.fileformat])
                           )
            sc.pl.umap(self.adata_subset, color=list(self.adata_subset.var_names), show=False, layer="scaled01",
                       legend_fontoutline=1, na_in_legend=False, s=30,
                       title=_, cmap='turbo', groups=[_],
                       save=".".join(["".join([str(self.tool), "_all"]), self.fileformat])
                       )
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    sc.pl.umap(self.adata_subset, color=_,
                               cmap=self.palette, legend_fontoutline=2, show=False, add_outline=False, frameon=False,
                               title="UMAP Plot",
                               s=50, save=".".join(["_".join([str(self.tool), _]), "pdf"]))

        elif self.runtime == 'clustering':
            self.createdir("/".join([self.output_folder, "".join(["Figures", self.tool])]))
            self.outfig = "/".join([self.output_folder, "".join(["Figures", self.tool])])

    def matrixplot(self):
        """
        Function for the generation of matrixplot sc.pl.matrixplot
        return:
        """
        if self.runtime != 'umap':
            sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                             dendrogram=True, vmin=-2, vmax=2, cmap='RdBu_r', layer="scaled",
                             show=False, swap_axes=False, return_fig=False,
                             save=".".join(["matrixplot_mean_z_score", self.fileformat]))
            sc.pl.matrixplot(self.adata_subset, list(self.adata_subset.var_names), "pheno_leiden",
                             dendrogram=True, cmap='Blues', standard_scale='var',
                             colorbar_title='column scaled\nexpression', layer="scaled",
                             swap_axes=False, return_fig=False,
                             show=False,
                             save=".".join(["matrixplot_column_scaled_expression", self.fileformat]))
        else:
            pass

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
        if self.runtime != 'clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        self.adata_subset.obs['pheno_leiden'], graph, Q = pg.cluster(self.adata_subset.X, k=int(self.k_coef),
                                                                     seed=42,
                                                                     clustering_algo="leiden",
                                                                     directed=True, primary_metric="euclidean",
                                                                     q_tol=0.05,
                                                                     prune=False, min_cluster_size=1,
                                                                     n_jobs=int(self.thread))
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden'].values
        self.adata.obs['Phenograph_cluster'] = self.adata_subset.obs['pheno_leiden'].values
        if self.runtime != 'clustering':
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        self.plot_frequency()
        self.matrixplot()
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
        if self.runtime != 'clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        p = parc.PARC(self.adata_subset.X, random_seed=42, knn=int(self.k_coef),
                      jac_std_global='median', jac_weighted_edges=False, distance='l2',
                      small_pop=10,
                      num_threads=int(self.thread))
        p.run_PARC()
        self.adata_subset.obs['pheno_leiden'] = [str(i) for i in p.labels]
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype(int) + 1
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['cluster'] = self.adata_subset.obs['pheno_leiden'].values
        self.adata.obs['Parc_cluster'] = self.adata_subset.obs['pheno_leiden'].values
        if self.runtime != 'clustering':
            self.embedding = self.runumap()
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
        self.generation_concatenate()
        self.plot_umap()
        self.plot_frequency()
        self.matrixplot()
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
        if self.runtime != 'clustering':
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset = self.adata[:, self.markertoinclude].copy()
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        else:
            if self.scanorama is True:
                self.adata_subset = self.correct_scanorama()
            else:
                self.adata_subset = self.adata[:, self.markertoinclude].copy()
                self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                                 zero_center=True, copy=True).X
                self.adata_subset.X = self.adata_subset.layers['scaled']
        self.adata_subset.X = self.adata_subset.layers['raw_value']
        tmpdf = self.adata_subset.to_df()
        tmpdf['Time'] = 0
        tmpdf.to_csv(tmp.name, header=True, index=False)
        tt = flowsom(tmp.name, if_fcs=False, if_drop=True)
        sample_df = tt.df
        if tt.df.shape[0] > 50000:
            self.flowsomshape= 20
        else:
            self.flowsomshape = 50
        tt.som_mapping(self.flowsomshape, self.flowsomshape, tt.df.shape[1],
                       sigma=2.5,
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
        output_tf_df = tt.tf_df  # new column added: category
        output_tf_df.set_index(self.adata.obs.index, inplace=True)
        # clustering ordering
        a = output_tf_df['category'].value_counts()
        a = a.reset_index()
        b = a[['index']].to_dict()
        # revert dict
        output_tf_df['category'] = output_tf_df['category'].map({v: k for k, v in b["index"].items()})
        # assign
        self.adata_subset.obs['pheno_leiden'] = output_tf_df['category'].values
        self.adata.obs['cluster'] = output_tf_df['category'].values
        # convert to category
        self.adata_subset.obs['pheno_leiden'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        self.adata.obs['Cluster_Flowsom'] = self.adata_subset.obs['pheno_leiden'].astype('category')
        #
        self.adata_subset.X = self.adata_subset.layers['scaled']
        self.embedding = self.runumap()
        if self.runtime != 'clustering':
            self.adata.obsm['X_umap'] = self.embedding
            self.adata_subset.obsm['X_umap'] = self.embedding
            for _ in list(self.adata_subset.obs['pheno_leiden'].unique()):
                if self.adata_subset[self.adata_subset.obs['pheno_leiden'].isin([_]), :].shape[0] < 20:
                    self.adata_subset = self.adata_subset[~self.adata_subset.obs['pheno_leiden'].isin([_]), :]
                    self.adata = self.adata[~self.adata.obs['Cluster_Flowsom'].isin([_]), :]
        self.generation_concatenate()
        self.plot_umap()
        self.plot_frequency()
        self.matrixplot()
        return self.adata

    def generation_concatenate(self):
        """

        """
        if self.runtime == 'full':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns=self.adata.var_names,
                                                index=self.adata.obs.index).astype(int),
                                   pd.DataFrame(self.adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'],
                                                index=self.adata.obs.index),
                                   right_index=True,
                                   left_index=True)
            pd.merge(self.tmp_df, self.adata.obs[['cluster',
                                                  'Sample', 'Cell_type',
                                                  'EXP',
                                                  'ID', 'Time_point',
                                                  'Condition']], left_index=True,
                     right_index=True).to_csv("/".join([self.output_folder,
                                                        "_ConcatenatedCells.".join(["_".join([self.analysis_name]),
                                                                                    "csv"])]),
                                              header=True, index=False)
        elif self.runtime == 'umap':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns=self.adata.var_names,
                                                index=self.adata.obs.index).astype(int),
                                   pd.DataFrame(self.embedding, columns=['UMAP_1', 'UMAP_2'],
                                                index=self.adata.obs.index),
                                   right_index=True,
                                   left_index=True)
            pd.merge(self.tmp_df, self.adata.obs[['Sample', 'Cell_type',
                                                  'EXP',
                                                  'ID', 'Time_point',
                                                  'Condition']], left_index=True,
                     right_index=True).to_csv("/".join([self.output_folder,
                                                        "_ConcatenatedCells.".join(["_".join([self.analysis_name]),
                                                                                    "csv"])]),
                                              header=True, index=False)
        elif self.runtime == 'clustering':
            self.tmp_df = pd.merge(pd.DataFrame(self.adata.X,
                                                columns=self.adata.var_names,
                                                index=self.adata.obs.index).astype(int),
                                   self.adata.obs[['cluster', 'Sample', 'Cell_type',
                                                   'EXP',
                                                   'ID', 'Time_point',
                                                   'Condition']],
                                   right_index=True,
                                   left_index=True)
            self.tmp_df.to_csv("/".join([self.output_folder,
                                         "_ConcatenatedCells.".join(["_".join([self.analysis_name]),
                                                                     "csv"])]),
                               header=True, index=False)

    def plot_frequency(self):
        """
        Plot barplot with frequency
        """
        if self.runtime != 'umap':
            fig, (ax1) = plt.subplots(1, 1, figsize=(17 / 2.54, 17 / 2.54))
            ax1 = self.adata_subset.obs.groupby("pheno_leiden")["Sample"].value_counts(
                normalize=True).unstack().plot.barh(
                stacked=True,
                legend=False,
                ax=ax1,
                color=self.palette)
            ax1.set_xlabel("Percentage Frequency")
            ax1.set_ylabel("Cluster")
            ax1.grid(False)
            ax1.legend(bbox_to_anchor=(1.2, 1.0))
            fig.savefig("/".join([self.outfig, ".".join(["ClusterFrequencyNormalized", self.fileformat])]),
                        dpi=self.dpi, bbox_inches='tight',
                        format=self.fileformat)
            #
            for _ in ['Sample', 'Cell_type', 'EXP', 'ID', 'Time_point', 'Condition']:
                if len(self.adata_subset.obs[_].unique()) > 1:
                    fig, (ax3) = plt.subplots(1, 1, figsize=(17 / 2.54, 17 / 2.54))
                    ax3 = self.adata_subset.T.var.groupby(_)["pheno_leiden"].value_counts(
                        normalize=True).unstack().plot.barh(
                        stacked=True,
                        legend=False,
                        color=self.palette,
                        ax=ax3,
                        fontsize=5)
                    ax3.set_xlabel("Cluster Percentage Frequency")
                    ax3.set_ylabel(_)
                    ax3.grid(False)
                    ax3.legend(bbox_to_anchor=(1.2, 1.0))
                    fig.savefig("/".join([self.outfig, ".".join(["".join([_, "ClusterFrequencyNormalized"]),
                                                                 self.fileformat])]),
                                dpi=self.dpi, bbox_inches='tight',
                                format=self.fileformat)

            #
            fig, (ax2) = plt.subplots(1, 1, figsize=(17 / 2.54, 17 / 2.54))
            ax2 = self.adata_subset.obs.groupby("pheno_leiden")["Sample"].value_counts(
                normalize=False).unstack().plot.barh(
                stacked=True,
                legend=False,
                ax=ax2,
                color=self.palette)
            ax2.set_xlabel("Relative Frequency")
            ax2.set_ylabel("Cluster")
            ax2.grid(False)
            ax2.legend(bbox_to_anchor=(1.2, 1.0))
            fig.savefig("/".join([self.outfig, ".".join(["ClusterFrequencyNotNormalized", self.fileformat])]),
                        dpi=self.dpi, bbox_inches='tight',
                        format=self.fileformat)
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
            self.tmp = self.adata[self.adata.obs['cluster'].isin([_])].to_df()
            self.tmp = self.tmp.astype(int)
            if self.runtime == 'full':
                self.tmp['UMAP_1'] = self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:, 0]
                self.tmp['UMAP_2'] = self.adata[self.adata.obs['cluster'].isin([_])].obsm['X_umap'][:, 1]
            if (self.tool == "Phenograph"):
                self.tmp['Phenograph'] = _
            elif (self.tool == "Parc"):
                self.tmp['Parc'] = _
            else:
                self.tmp['Flowsom'] = _
            self.tmp.to_csv("/".join([self.output_folder, "".join(["FCScluster",
                                                                   self.tool]), "".join([self.analysis_name, "_",
                                                                                         str(_), ".csv"])]),
                            header=True, index=False)

    def groupbysample(self):
        """
        Function for generation of csv with different clusters
        """
        # make dir
        self.createdir("/".join([self.output_folder, "".join(["FCSsample", self.tool])]))
        # tmp df
        self.tmp = self.adata.to_df()
        self.tmp = self.tmp.astype(int)
        # change index
        self.tmp.set_index(self.adata.obs['Sample'], inplace=True)
        if self.runtime != 'umap':
            self.adata.obs['cluster'] = self.adata.obs['cluster'].astype(int)
            self.createdir("/".join([self.output_folder, "".join(["Cluster_frequency", self.tool])]))
        # create columns
        if self.runtime != 'clustering':
            self.tmp['UMAP_1'] = self.adata.obsm['X_umap'][:, 0]
            self.tmp['UMAP_2'] = self.adata.obsm['X_umap'][:, 1]
        if self.runtime != 'umap':
            self.tmp[self.tool] = self.adata.obs['cluster'].values
            self.tmp["cluster"] = self.adata.obs['cluster'].values
            # get unique filenames
            unique_filename = self.adata.obs['Sample'].unique()
            # get unique number of cluster
            unique_Phenograph = self.adata.obs['cluster'].unique()
            #
            dfCounts = pd.DataFrame(index=range(min(unique_Phenograph), max(unique_Phenograph) + 1))
            # generate Tot_percentage file
            for i in range(len(unique_filename)):
                dfCounts[unique_filename[i]] = self.tmp.loc[self.tmp.index == unique_filename[i]].cluster.value_counts(
                    normalize=True).reindex(self.tmp.cluster.unique(), fill_value=0)
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
                dfCounts[unique_filename[i]] = self.tmp.loc[self.tmp.index == unique_filename[i]].cluster.value_counts().reindex(
                    self.tmp.cluster.unique(), fill_value=0)
            # save
            dfCounts.index.name = 'Cluster'
            dfCounts.to_csv("/".join(["/".join([self.output_folder, "".join(["Cluster_frequency", self.tool])]),
                                      "".join(["Tot_counts", ".csv"])]))
            del self.tmp['cluster']
            # save samples
            for i in range(len(unique_filename)):
                dfCounts[unique_filename[i]] = self.tmp.loc[self.tmp.index == unique_filename[i]].to_csv(
                    "/".join([self.output_folder, "".join(["FCSsample", self.tool]),
                              "".join([str(unique_filename[i]), "_", self.analysis_name,
                                       ".csv"])]),
                    header=True, index=False)
        else:
            unique_filename = self.adata.obs['Sample'].unique()
            for i in range(len(unique_filename)):
                self.tmp.loc[self.tmp.index == unique_filename[i]].to_csv(
                    "/".join([self.output_folder, "".join(["FCSsample", self.tool]),
                              "".join([str(unique_filename[i]), "_", self.analysis_name,
                                       ".csv"])]),
                    header=True, index=False)

    def runtimeumap(self):
        """
        Function for execution of phenograph analysis
        :return:
        """
        self.log.info("Part2: UMAP dimensional reduction")
        # self.log.info("Markers used for Phenograph clustering:")
        self.adata_subset = self.adata[:, self.markertoinclude].copy()
        self.log.info(self.adata_subset.var_names)
        self.log.info("Markers excluded for UMAP computation:")
        self.log.info(self.marker_array)
        self.adata_subset.layers['scaled'] = sc.pp.scale(self.adata_subset, max_value=6,
                                                         zero_center=True, copy=True).X
        #self.adata_subset.X = self.adata_subset.layers['scaled']
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
        self.log.info("Part4: Output Generation")
        self.scaler = MinMaxScaler(feature_range=(0, 1))
        if self.runtime != 'umap':
            if self.tool == "Phenograph" or self.tool == "Parc":
                self.adata.obs[self.tool + "_" + str(self.k_coef)] = self.adata.obs['cluster'].astype("str")
                del self.adata.obs['cluster']
                del self.adata.obs[self.tool + "_" + str(self.k_coef)]
                self.adata.obs["".join([str(self.tool), "_cluster"])] = self.adata.obs[
                    "".join([str(self.tool), "_cluster"])].astype('category')
                self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
                self.adata.X = self.adata.layers['scaled01']
                self.adata.layers['scaled01'] =scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
                self.adata.write(
                    "/".join([self.output_folder, ".".join(["_".join([self.analysis_name, self.k_coef]), "h5ad"])]))
            elif self.tool == "Flowsom":
                del self.adata.obs['cluster']
                self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
                self.adata.X = self.adata.layers['scaled01']
                self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
                self.adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "h5ad"])]))
        else:
            self.adata.layers['scaled01'] = self.scaler.fit_transform(self.adata.layers['raw_value'])
            self.adata.X = self.adata.layers['scaled01']
            self.adata.layers['scaled01'] = scipy.sparse.csr_matrix(self.adata.layers['scaled01'])
            self.adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name]), "h5ad"])]))
