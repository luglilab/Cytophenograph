from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

print(__doc__)


data = pd.read_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/Jasper/CD4_Ag_Phenograph/30/Exported_FJ10_CD4_Ag_k30_concatenated.txt",
                   sep="\t", header=0)

print(data.head())

data = concdf.copy()
for i in range(len(marker)):
    data.drop(marker[i], axis=1, inplace=True)
x = data.values




#
# y = data[["Phenograph"]].values.flatten()
# X = data.drop(columns=['Phenograph']).values
# #print(X.shape)
#
# range_n_clusters = [23,11]
#
# for n_clusters in range_n_clusters:
#     # Create a subplot with 1 row and 2 columns
#     fig, (ax1, ax2) = plt.subplots(1, 2)
#     fig.set_size_inches(18, 7)
#
#     # The 1st subplot is the silhouette plot
#     # The silhouette coefficient can range from -1, 1 but in this example all
#     # lie within [-0.1, 1]
#     ax1.set_xlim([-0.1, 1])
#     # The (n_clusters+1)*10 is for inserting blank space between silhouette
#     # plots of individual clusters, to demarcate them clearly.
#     ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
#
#     # Initialize the clusterer with n_clusters value and a random generator
#     # seed of 10 for reproducibility.
#     cluster_labels = y
#     silhouette_avg = silhouette_score(X, cluster_labels)
#     # clusterer = KMeans(n_clusters=n_clusters, random_state=10)
#     # cluster_labels = clusterer.fit_predict(X)
#     #
#     # # The silhouette_score gives the average value for all the samples.
#     # # This gives a perspective into the density and separation of the formed
#     # # clusters
#     # silhouette_avg = silhouette_score(X, cluster_labels)
#     print("For n_clusters =", n_clusters,
#           "The average silhouette_score is :", silhouette_avg)
#
#     # Compute the silhouette scores for each sample
#     sample_silhouette_values = silhouette_samples(X, cluster_labels)
#
#     y_lower = 10
#     for i in range(n_clusters):
#         # Aggregate the silhouette scores for samples belonging to
#         # cluster i, and sort them
#         ith_cluster_silhouette_values = \
#             sample_silhouette_values[cluster_labels == i]
#
#         ith_cluster_silhouette_values.sort()
#
#         size_cluster_i = ith_cluster_silhouette_values.shape[0]
#         y_upper = y_lower + size_cluster_i
#
#         color = cm.nipy_spectral(float(i) / n_clusters)
#         ax1.fill_betweenx(np.arange(y_lower, y_upper),
#                           0, ith_cluster_silhouette_values,
#                           facecolor=color, edgecolor=color, alpha=0.7)
#
#         # Label the silhouette plots with their cluster numbers at the middle
#         ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
#
#         # Compute the new y_lower for next plot
#         y_lower = y_upper + 10  # 10 for the 0 samples
#
#     ax1.set_title("The silhouette plot for the various clusters.")
#     ax1.set_xlabel("The silhouette coefficient values")
#     ax1.set_ylabel("Cluster label")
#
#     # The vertical line for average silhouette score of all the values
#     ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
#
#     ax1.set_yticks([])  # Clear the yaxis labels / ticks
#     ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
#
#     # 2nd Plot showing the actual clusters formed
#     colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
#     ax2.scatter(x[:,0], x[:,1], lw=0, s=40,alpha=0.3, c=palette[colors.astype(np.int)])
#     ax2.xlim(-25, 25)
#     ax2.ylim(-25, 25)
#     ax2.axis('off')
#     ax2.axis('tight')
#     # add the labels for each digit corresponding to the label
#     txts = []
#     for i in range(num_classes):
#
#         # Position of each label at median of data points.
#
#         xtext, ytext = np.median(x[colors == i, :], axis=0)
#         txt = ax2.text(xtext, ytext, str(i), fontsize=24)
#         txt.set_path_effects([
#             PathEffects.Stroke(linewidth=5, foreground="w"),
#             PathEffects.Normal()])
#         txts.append(txt)
#
#
#
#     # ax2.set_title("The visualization of the clustered data.")
#     # ax2.set_xlabel("Feature space for the 1st feature")
#     # ax2.set_ylabel("Feature space for the 2nd feature")
#
#     plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
#                   "with n_clusters = %d" % n_clusters),
#                  fontsize=14, fontweight='bold')
#
# plt.show()