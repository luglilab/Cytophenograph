import pytest
import numpy as np
from scipy import sparse
from sklearn import datasets

class ClusterDataset:
    """
    Simple class to return dataset information and share between test cases
    """

    def __init__(self):
        """
        Object constructor

        Returns:
            ClusterDataset: An example dataset for use with PhenoGraph testing
        """

        self.dense_data = np.array(
            [[  8.29188122,   6.77467787],
            [  8.11595558,   6.29969468],
            [  8.41932782,   7.18708683],
            [  7.87751289,   7.05695946],
            [  8.15498076,   7.2317535 ],
            [  8.57681   ,   7.30907322],
            [  7.84036909,   7.04545501],
            [ -9.22407124, -10.09577064],
            [ -8.97447722, -10.29175587],
            [ -9.47952697, -10.26263977],
            [ -9.13671857,  -9.98270161],
            [ -9.17054197, -10.26074704],
            [ -9.20421828, -10.40169291],
            [ -9.28332716, -10.29473222],
            [ -5.9341068 , -19.0755188 ],
            [ -5.83337714, -19.25045898],
            [ -6.12433156, -19.51021562],
            [ -5.60446077, -19.16554577],
            [ -5.56675018, -19.50978007],
            [ -5.93608833, -18.91981749],
            [ -5.7953458 , -19.04092524],
            [  8.00608704,   7.40402037],
            [  8.12470637,   7.1818027 ],
            [  8.36344614,   7.17558329],
            [ -9.19624784, -10.08752065],
            [ -9.41859495, -10.31316491],
            [ -9.30012735, -10.32240619],
            [ -6.20397103, -19.56631818],
            [ -5.83048038, -19.3580962 ],
            [ -6.21581103, -19.22808712]])
        self.sparse_data = sparse.coo_matrix(self.dense_data)
        self.initial_labels = np.array(
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])
        self.actual_labels = np.array(
            [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2])
        self.train_set = [self.dense_data[i:i+3] for i in np.arange(21, 30, 3)]
        self.test_set = self.dense_data[:21]
        self.kNN_matrix = np.array(
            [[23,  2, 22,  4,  3],
            [ 0,  3,  6, 22, 23],
            [23,  5,  4, 22,  0],
            [ 6, 22,  4, 21, 23],
            [22, 23, 21,  2,  3],
            [ 2, 23,  4, 22, 21],
            [ 3, 22,  4, 21,  0],
            [24, 10, 11, 13, 26],
            [11, 12, 24, 13,  7],
            [25, 26, 13,  7, 12],
            [24,  7, 11, 13,  8],
            [13, 26, 12,  7, 24],
            [26, 13, 11, 25,  8],
            [26, 11, 12, 25,  9],
            [20, 19, 15, 28, 29],
            [28, 14, 20, 17, 19],
            [27, 29, 28, 15, 14],
            [20, 15, 28, 14, 18],
            [28, 17, 15, 20, 16],
            [14, 20, 15, 17, 29],
            [14, 19, 15, 17, 28],
            [ 4, 22,  3,  6, 23],
            [ 4, 23, 21,  3,  2],
            [ 2,  4, 22,  5,  0],
            [ 7, 10, 11, 13, 26],
            [ 9, 26, 13, 12, 11],
            [13, 25, 12, 11,  9],
            [16, 29, 28, 15, 14],
            [15, 17, 14, 18, 20],
            [16, 14, 27, 15, 28]])
        self.sparse_affinity_row = np.array(
            [ 0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,
            3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  6,  6,  6,  6,
            6,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9, 10,
            10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13,
            13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16,
            17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20,
            20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23,
            23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 27,
            27, 27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29],
            )
        self.sparse_affinity_col = np.array(
            [23,  2, 22,  4,  3,  0,  3,  6, 22, 23, 23,  5,  4, 22,  0,  6, 22,
            4, 21, 23, 22, 23, 21,  2,  3,  2, 23,  4, 22, 21,  3, 22,  4, 21,
            0, 24, 10, 11, 13, 26, 11, 12, 24, 13,  7, 25, 26, 13,  7, 12, 24,
            7, 11, 13,  8, 13, 26, 12,  7, 24, 26, 13, 11, 25,  8, 26, 11, 12,
            25,  9, 20, 19, 15, 28, 29, 28, 14, 20, 17, 19, 27, 29, 28, 15, 14,
            20, 15, 28, 14, 18, 28, 17, 15, 20, 16, 14, 20, 15, 17, 29, 14, 19,
            15, 17, 28,  4, 22,  3,  6, 23,  4, 23, 21,  3,  2,  2,  4, 22,  5,
            0,  7, 10, 11, 13, 26,  9, 26, 13, 12, 11, 13, 25, 12, 11,  9, 16,
            29, 28, 15, 14, 15, 17, 14, 18, 20, 16, 14, 27, 15, 28],
            )
        self.sparse_affinity_data = np.array(
            [0.42857143, 0.42857143, 0.66666667, 0.66666667, 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.25      , 0.25      ,
            0.66666667, 0.42857143, 0.25      , 0.25      , 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.66666667, 0.25      ,
            0.66666667, 0.25      , 0.42857143, 0.25      , 0.42857143,
            0.42857143, 0.42857143, 0.66666667, 0.66666667, 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.42857143, 0.42857143,
            0.66666667, 0.42857143, 0.42857143, 0.25      , 0.25      ,
            0.66666667, 0.25      , 0.42857143, 0.25      , 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.25      , 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.11111111, 0.66666667,
            0.25      , 0.25      , 0.25      , 0.42857143, 0.42857143,
            0.42857143, 0.42857143, 0.25      , 0.42857143, 0.25      ,
            0.66666667, 0.25      , 0.42857143, 0.66666667, 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.25      , 0.25      ,
            0.42857143, 0.42857143, 0.66666667, 0.42857143, 0.42857143,
            0.66666667, 0.66666667, 0.25      , 0.25      , 0.42857143,
            0.42857143, 0.42857143, 0.66666667, 0.42857143, 0.42857143,
            0.42857143, 0.42857143, 0.42857143, 0.42857143, 0.25      ,
            0.42857143, 0.42857143, 0.42857143, 0.42857143, 0.25      ,
            0.42857143, 0.42857143, 0.66666667, 0.42857143, 0.42857143,
            0.42857143, 0.42857143, 0.66666667, 0.42857143, 0.25      ,
            0.66666667, 0.25      , 0.42857143, 0.42857143, 0.25      ,
            0.66666667, 0.25      , 0.25      , 0.42857143, 0.42857143,
            0.66666667, 0.42857143, 0.42857143, 0.25      , 0.25      ,
            0.42857143, 0.66666667, 0.66666667, 0.42857143, 0.42857143,
            0.66666667, 0.66666667, 0.42857143, 0.25      , 0.42857143,
            0.66666667, 0.66666667, 0.25      , 0.25      , 0.42857143,
            0.42857143, 0.66666667, 0.25      , 0.42857143, 0.42857143,
            0.66666667, 0.25      , 0.66666667, 0.25      , 0.25      ])
        self.sparse_affinity_matrix = sparse.csr_matrix((self.sparse_affinity_data, (self.sparse_affinity_row, self.sparse_affinity_col)))
        self.random_walk_P = np.array(
            [[1., 0., 0.],
            [1., 0., 0.],
            [1., 0., 0.],
            [1., 0., 0.],
            [1., 0., 0.],
            [1., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.]])
        self.coo_affinity_matrix = sparse.coo_matrix(self.sparse_affinity_matrix)
        self.communities = np.array(
            [0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 1])
        self.Q_louvain = 0.666593
        self.Q_leiden = None


class GeneratedClusterDataset(ClusterDataset):
    
    """
    Simple class wrapping Scikit Learn's classification function
    Intended to be used for testing Phenograph with pytest fixtures
    
    Attributes:
        n_clusters (int): No. of distinct clusters
        n_points (int): No. of points to include in training set
        n_train_pts (int): No. of points to include in training set
        n_dims (int, optional): No. of dimensions of dataset. Defaults to 2.
        test_set (numpy array): Points included in test set
        train_set (list of numpy array): Points included in training set, ordered by label
        test_labels (list of int): Correct labels for points in test set
    """

    def __init__(self, n_clusters:int, n_points:int, n_train_pts:int, n_dims:int = 2):
        """
        Object constructor

        Args:
            n_clusters (int): No. of distinct clusters
            n_points (int): No. of points to generate per cluster
            n_train_pts (int): No. of points to include in training set
            n_dims (int, optional): No. of dimensions of dataset. Defaults to 2.
        """            

        assert(n_points > n_train_pts)
        assert(n_dims >= 2)

        self.n_clusters = n_clusters
        self.n_points = n_points
        self.n_train_pts = n_train_pts
        self.n_dims = n_dims

        data, clusters = datasets.make_blobs(
            self.n_points*self.n_clusters,
            self.n_dims,
            centers = self.n_clusters,
            cluster_std=0.25,
            center_box=(-20, 20)
            )

        self.train_set = []
        test_indices = set(range(self.n_points*self.n_clusters))
        for n in range(self.n_clusters):
            train_indices = np.argwhere(clusters == n).flatten()[:self.n_train_pts]
            self.train_set.append(data[train_indices])
            test_indices -= set(train_indices)

        test_indices = list(test_indices)
        self.test_set = data[test_indices]
        self.test_labels = clusters[test_indices]   

@pytest.fixture(scope="module")
def cluster_dataset_generated():
    cluster_test = GeneratedClusterDataset(3, 25, 24, 2)
    return cluster_test

@pytest.fixture(scope="module")
def cluster_dataset_fixed():
    cluster_test = ClusterDataset()
    return cluster_test