from phenograph.cluster import run_louvain, run_leiden
from scipy.sparse import coo_matrix
from numpy.testing import assert_array_equal
import numpy as np


def compare_communities(A, B):
    communities_A, idx = np.unique(A, return_index=True)
    communities_B, idx = np.unique(B, return_index=True)
    if not np.array_equal(communities_A, communities_B):
        return False
    else:
        communities_A = A[np.sort(idx)]
        communities_B = B[np.sort(idx)]
        offset = communities_A.max() + 1
        for i, comm in enumerate(communities_A):
            A = np.where(A == comm, i + offset, A)
        for j, comm in enumerate(communities_B):
            B = np.where(B == comm, j + offset, B)
        return np.array_equal(A, B)


def test_run_louvain(cluster_dataset_fixed):
    q_tol = 1e-3
    louvain_time_limit = 2000
    graph = cluster_dataset_fixed.coo_affinity_matrix
    communities, Q = run_louvain(graph, q_tol, louvain_time_limit)
    assert compare_communities(communities, cluster_dataset_fixed.communities)
    assert Q == cluster_dataset_fixed.Q_louvain


def test_run_leiden(cluster_dataset_fixed):
    directed = True
    partition_type = None
    resolution_parameter = 1
    n_iterations = -1
    seed = None
    use_weights = True
    kargs = dict()
    graph = cluster_dataset_fixed.coo_affinity_matrix
    communities, Q = run_leiden(
        graph,
        directed,
        partition_type,
        resolution_parameter,
        n_iterations,
        seed,
        use_weights,
        kargs,
    )
    assert compare_communities(communities, cluster_dataset_fixed.communities)
    assert Q == cluster_dataset_fixed.Q_leiden
