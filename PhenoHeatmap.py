from numpy import nanmin
import matplotlib.pyplot as plt
plt.switch_backend('agg') # prevent show plot bug
from seaborn import heatmap, clustermap
import pandas as pd
import seaborn as sns
import numpy as np
sns.set()



def plot_heatmap(frame, save_to_path, y_axis):

	land = frame.shape[1] > frame.shape[0]
	plt.figure(figsize=get_A4_dimens(is_landscape=land))
	square = True if frame.shape[0] < 25 else False # ensure fits A4

	hm = heatmap(frame, cmap='jet', xticklabels=True, yticklabels=True, square=square).get_figure()
	plt.xlabel('marker')
	plt.ylabel(y_axis)
	plt.xticks(rotation=90)
	plt.yticks(rotation=0, fontsize=6)

	save_pdf(save_to_path, hm)
	print('Export {}'.format(save_to_path), flush=True)
	plt.close('all')


df1 = pd.read_csv("/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/outfoldertest/testGiorgiaFCS_cluster_info.csv",sep="\t",header=0)
df1 = df1.drop(columns=['Unnamed: 0','cluster'])


test = sns.heatmap(df1.T,
                 vmin=0,
                 vmax=1.0,
                 linewidths=0.004,
                 linecolor="black",
                 cmap='Spectral_r')
test = test.get_figure()
test.savefig("heatmap.png")

