from phenograph.core import parallel_jaccard_kernel, neighbor_graph
from numpy.testing import assert_array_equal, assert_array_almost_equal


def test_parallel_jaccard_kernel(cluster_dataset_fixed):
    actual = cluster_dataset_fixed
    row, col, data = parallel_jaccard_kernel(actual.kNN_matrix)
    assert_array_equal(row, actual.sparse_affinity_row)
    assert_array_equal(col, actual.sparse_affinity_col)
    assert_array_almost_equal(data, actual.sparse_affinity_data)


def test_neighbor_graph(cluster_dataset_fixed):
    assert True