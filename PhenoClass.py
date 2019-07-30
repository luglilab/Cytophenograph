import glob
import os
import sys
import pandas as pd
import phenograph as pg
from openTSNE import TSNE
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.cm as cm
import numpy as np
from sklearn.metrics import silhouette_samples, silhouette_score
import random
import subprocess
from numpy import unique
import umap
import warnings
warnings.filterwarnings('ignore')


class cytopheno(object):

    def __init__(self, optparseinstance):




