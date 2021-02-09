import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import anndata
import glob
import os
import sys
import numpy as np
import pandas as pd
import re
# import seaborn as sb
import matplotlib.pyplot as plt
#import scanpy.external as sce
import phenograph as pg
import scanpy as sc
#from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm
import parc
import umap
import seaborn as sns
import matplotlib.patheffects as PathEffects
#from contextlib import contextmanager,redirect_stderr,redirect_stdout
# from os import devnull
import logging
#warnings.filterwarnings('ignore')
#warnings.filterwarnings("ignore", message="Numerical issues were encountered ")


class Cytophenograph:
    def __init__(self, info_file, input_folder, output_folder, k_coef, marker_list, analysis_name, thread, tsne,tool):
        self.info_file = info_file
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.k_coef = k_coef
        self.marker_list = marker_list
        self.analysis_name = analysis_name
        self.thread = thread
        self.tsne = tsne
        self.tool = tool
        sys.stdout = open("/".join([self.output_folder,"log.txt"]), 'a')
        self.log = logging.getLogger()
        self.log.setLevel(logging.DEBUG)
        format = logging.Formatter("%(asctime)s %(threadName)-11s %(levelname)-10s %(message)s")
        #
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        self.log.addHandler(ch)
        # 
        fh = logging.FileHandler("/".join([self.output_folder,"log.txt"]), "w")
        fh.setFormatter(format)
        self.log.addHandler(fh)

        #self.log = logging.getLogger()
        #self.log.setLevel(logging.INFO)
        #self.log.basicConfig(stream=sys.stdout, level=logging.DEBUG)
        #self.filehandler = logging.StreamHandler(sys.stdout) 
        #self.filehandler = logging.FileHandler("/".join([self.output_folder,"log.txt"]), "w")
        #self.stdout_handler = logging.StreamHandler(sys.stdout)
        #self.filehandler.setLevel(logging.DEBUG)
        #formatter = logging.Formatter( "%(asctime)s %(threadName)-11s %(levelname)-10s %(message)s")
        #self.filehandler.setFormatter(formatter)
        #self.stdout_handler.setFormatter(formatter)
        #self.log.addHandler(self.filehandler)
        #self.log.addHandler(self.stdout_handler)
        self.log.info("Name of this analysis: {}".format(marker_list))
        self.log.info("Input folder: {}".format(input_folder))
        self.log.info("Output folder: {}".format(output_folder))
        self.log.info("Info file: {}".format(info_file))
        self.log.info("Phenograph K-coef: {}".format(k_coef))
        self.log.info("Marker list file: {}".format(marker_list))
        self.log.info("Clustering tool option: {}".format(tool))
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass  

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
                                ann_tmp.obs['Cell_type'] = cell_type
                                #
                                exp = info_file['EXP'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['EXP'] = exp
                                #
                                id = info_file['ID'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['ID'] = id
                                #
                                time_point = info_file['Time_point'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Time_point'] = time_point
                                #
                                condition = info_file['Condition'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Condition'] = condition
                                #
                                count = info_file['Count'].loc[info_file['Sample']== pandas_df_list[i].index[0][:-2]]
                                ann_tmp.obs['Count'] = count
                                anndata_list.append(ann_tmp)
                            else:
                                print("Error, this file {0} is not in the column Sample of Infofile. \n Please check sample name and Infofile".format(pandas_df_list[i].index[0][:-2]))
                                # print(pandas_df_list[i].index[0][:-2])
                                # print(Sample_list)
                                sys.exit(1)
                        tmp = anndata_list[0]
                        anndata_list.pop(0)
                        adata_conc = tmp.concatenate(anndata_list)
                    except (ValueError, Exception):
                        print("Error. Please check Info File Header or CSV header.")
                        self.log.error("Error. Please check Info File Header or CSV header.")   
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
        communities, graph, Q = pg.cluster(data.values, k=int(self.k_coef),directed=False,
        prune=False,min_cluster_size=1,n_jobs=int(self.thread))
        # create dataframe with Phenograph output
        dfPheno = pd.DataFrame(communities)
        # shift of one unit the name of cluster
        dfPheno["Phenograph"] = dfPheno[0] + 1
        # remove first column
        dfPheno = dfPheno.drop(columns=[0], axis=1)
        dfPheno.set_index(adata.obs.index, inplace=True)
        adata.obs['cluster'] = dfPheno
        adata.obs['Phenograph_cluster'] = dfPheno
        #print(adata.obs['Sample'].unique())
        reducer = umap.UMAP(random_state=42, n_neighbors=10, min_dist=0.001)
        embedding = reducer.fit_transform(data.values)
        adata.obsm['X_umap'] = embedding
        if self.tsne:
            from openTSNE import TSNE
            tsne = TSNE(n_components=2, perplexity=30, learning_rate=200,
                        n_jobs=self.thread, initialization="pca", metric="euclidean",
                        early_exaggeration_iter=250, early_exaggeration=12, n_iter=1000,
                        neighbors="exact", negative_gradient_method="bh")
            embedding = tsne.fit(data.values)
            adata.obsm['X_tsne'] = embedding
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
        #print(adata.obs['Sample'].unique())
        if self.tsne:
            tsne = TSNE(n_components=2, perplexity=30, learning_rate=200,
                        n_jobs=self.thread, initialization="pca", metric="euclidean",
                        early_exaggeration_iter=250, early_exaggeration=12, n_iter=1000,
                        neighbors="exact", negative_gradient_method="bh")
            embedding = tsne.fit(data.values)
            adata.obsm['X_tsne'] = embedding
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
            tmp['Phenograph'] = _
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
        # save samples
        for i in range(len(unique_filename)):
            dfCounts[unique_filename[i]] = tmp.loc[tmp.index == unique_filename[i]].to_csv(
                "/".join([self.output_folder, "".join(["FCSsample",tool]), "".join([str(unique_filename[i]), "_", self.analysis_name,
                                                                    ".csv"])]),
                header=True, index=False)


    def suppress_stdout_stderr(self):
        """A context manager that redirects stdout and stderr to devnull"""
        with open(devnull, 'w') as fnull:
            with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
                yield (err, out)

    def tsne_umap_plot(self,x, cluster, kind):
        """
        Tsne_1  Tsne_2
        :param x:
        :param cluster:
        :param kind:
        :return:
        """
        colors = cluster
        # x = x.values
        # choose a color palette with seaborn.
        num_classes = len(np.unique(colors))
        palette = np.array(sns.color_palette("hls", num_classes + 1))  # outbound error
        # create a scatter plot.
        f = plt.figure(figsize=(8, 8))
        ax = plt.subplot(aspect='equal')
        scatter = ax.scatter(x[:, 0], x[:, 1], lw=0, s=15, alpha=0.5, c=palette[colors])  # flatten np array
        plt.xlim(-25, 25)
        plt.ylim(-25, 25)
        ax.axis('off')
        ax.axis('tight')
        # add the labels for each digit corresponding to the label
        txts = []
        for i in range(1, num_classes + 1):
            # Position of each label at median of data points.

            xtext, ytext = np.median(x[colors == i, :], axis=0)
            txt = ax.text(xtext, ytext, str(i), fontsize=24)
            txt.set_path_effects([
                PathEffects.Stroke(linewidth=5, foreground="w"),
                PathEffects.Normal()])
            txts.append(txt)

        # produce a legend with the unique colors from the scatter
        plt.title("Data embedded into two dimensions by {}".format(kind), fontsize=12)
        plt.savefig("/".join([self.output_folder, "".join([self.analysis_name, "_", kind, ".pdf"])]), format="pdf")
        #ax.set_rasterized(True)
        #plt.savefig("/".join([self.output_folder, "".join([self.analysis_name, "_", kind, ".eps"])]), format="eps")
        # plt.show()
        # return f, ax, scatter, txts
        return palette

    def exporting(self, adata):
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
        adata.obs["Phenograph_"+str(self.k_coef)] = adata.obs['cluster'].astype("str")
        del adata.obs['cluster']
        adata.write("/".join([self.output_folder, ".".join(["_".join([self.analysis_name, self.k_coef]), "h5ad"])]))






















def validationplot(markertoexclude, adata,outfold,name):
    """
    :param marker:
    :param alldf:
    file_name,event,FSC-A,FSC-H,FSC-W,SSC-A,SSC-H,Comp-APC-A,......,Tbet,Time,,Tsne_1,Tsne_2,Phenograph
    :param outfold:
    :param name:
    :return:
    """
    marker = adata.var_names.to_list()
    markertoinclude = [i for i in marker if i not in markertoexclude]
    data = adata[:,markertoinclude].to_df()
    
    # convert series in numpy array ordered
    y = adata.obs["cluster"].values 
    # convert from matrix to values
    X = data.values
    #
    n_clusters =len(np.unique(adata.obs['cluster']))
    #
    f = plt.figure(figsize=(18, 7))
    plt.style.use('seaborn-white')
    plt.xlim(-1, 1)
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    plt.ylim(0, len(X) + (n_clusters + 1) * 10)
    #
    cluster_labels = y
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(1,n_clusters+1):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters+1)
        plt.fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples


    plt.xlabel("The silhouette coefficient values ")
    plt.ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    plt.axvline(x=silhouette_avg, color="red", linestyle="--")

    plt.yticks([])  # Clear the yaxis labels / ticks
    plt.xticks([-1,-0.8,-0.6,-0.4,-0.2,-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    plt.title(("Silhouette analysis for Phenograph clustering on {0} data with n_clusters = {1}.\n"
               " The average silhouette_score is {2}".format(
                      name, n_clusters, silhouette_avg)),
              fontsize=14, fontweight='bold')
    plt.savefig("/".join([outfold, "".join([name, "_validation"])]))


def exporting(adata,outfold,name,kcoev):
    old_names = adata.var_names
    new_names = []
    for _ in range(len(old_names)):
        if old_names[_].startswith("Comp-"):
            new_names.append(old_names[_].split(":: ")[-1])
        else:
            new_names.append(old_names[_])
    adata.var = pd.DataFrame(old_names,new_names)
    del adata.var[0]
    adata.var['original_names'] = old_names
    adata.obs[str(kcoev)] = adata.obs['cluster'].astype("str")
    del adata.obs['cluster']
    adata.write("/".join([outfold,".".join(["_".join([name,kcoev]),"h5ad"])])) 

def lineplot(markertoexclude, adata,outfold,name,pal):
    marker = adata.var_names.to_list()
    markertoinclude = [i for i in marker if i not in markertoexclude]
    data = adata[:,markertoinclude].to_df()
    data["cluster"] = adata.obs["cluster"].values
    mmm = data.groupby(['cluster']).median()
    plt.figure(figsize=(15,10))
    x = [i for i in range(len(mmm.columns))]
    my_xticks = mmm.columns
    plt.xticks(x, my_xticks)
    for i in range(len(mmm)):
        plt.plot(x, mmm.iloc[i],c=pal[i+1])
    for i in range(1,len(mmm.columns)):
        plt.axvline(x=i, linestyle="dashdot", linewidth=0.5,c="lightgray")
    plt.xticks(rotation=90)
    plt.ylabel('Median channel value')
    plt.legend(mmm.index)
    plt.savefig("/".join([outfold, "".join([name, "_lineplot"])]))
