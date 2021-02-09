import multiprocessing as mp
import os
import re
import time
import uuid
from typing import Union, Optional, Type, Tuple

import igraph as ig
import leidenalg
import numpy as np
from leidenalg.VertexPartition import MutableVertexPartition
from scipy import sparse as sp
from scipy.sparse.base import spmatrix

from phenograph.core import (
    gaussian_kernel,
    parallel_jaccard_kernel,
    jaccard_kernel,
    find_neighbors,
    neighbor_graph,
    graph2binary,
    runlouvain,
)


def chunk_clusters(cl):
    for i in range(0, np.unique(cl).size, 5000):
        yield np.unique(cl)[i : i + 5000]


def yield_clusters(cl, ch):
    for i in ch:
        yield cl == i


def get_sizes(func, args):
    results = func(*args)
    return [np.count_nonzero(res) for res in results]


def sort_by_size(clusters: np.array, min_size: int = 10, n_jobs: int = -1) -> np.array:
    """\
    Relabel clustering in order of descending cluster size.
    New labels are consecutive integers beginning at 0.
    Clusters that are smaller than min_size are assigned to -1.

    Parameters
    ----------
    clusters
        Array of clusters to be sorted by size
    min_size
        Clusters smaller than this threshold are considered outliers and are assigned to
        -1 in the cluster labels
    n_jobs
        Number of concurrently running workers. If 1 is given, no parallelism is used.
        If set to -1, all CPUs are used. For n_jobs below -1, `n_cpus + 1 + n_jobs` are
        used.

    Returns
    -------
    Sorted array of clusters
    """
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    if n_jobs < -1:
        n_jobs = mp.cpu_count() + 1 + n_jobs
    p = mp.Pool(n_jobs)
    sizes = []
    ch_clust = chunk_clusters(clusters)
    TASKS = [(yield_clusters, (clusters, i)) for i in ch_clust]
    results = [p.apply_async(get_sizes, t) for t in TASKS]
    for res in results:
        sizes.extend(res.get())

    p.close()
    p.join()

    o = np.argsort(sizes)[::-1]
    my_dict = {c: i for i, c in enumerate(o) if sizes[c] > min_size}
    my_dict.update({c: -1 for i, c in enumerate(o) if sizes[c] <= min_size})

    relabeled = np.vectorize(my_dict.get)(clusters)

    return relabeled


def run_leiden(
    graph: sp.coo_matrix,
    directed: bool,
    partition_type: Optional[Type[MutableVertexPartition]],
    resolution_parameter: float,
    n_iterations: int,
    seed: Optional[int],
    use_weights: bool,
    kargs,
) -> Tuple[np.ndarray, float]:
    """
    Wrapper for leiden community detection

    Args:
        graph (sp.coo_matrix): Affinity matrix
        directed (bool): See below in 'cluster()'
        partition_type (Optional[Type[MutableVertexPartition]]): See below in 'cluster()'
        resolution_parameter (float): See below in 'cluster()'
        n_iterations (int): See below in 'cluster()'
        seed (Optional[int]): See below in 'cluster()'
        use_weights (bool): See below in 'cluster()'
        kargs: See below in 'cluster()'

    Returns:
        communities, Q (Tuple[np.ndarray, float]): See below in 'cluster()'
    """

    # convert resulting graph from scipy.sparse.coo.coo_matrix to Graph object
    # get indices of vertices
    edgelist = np.vstack(graph.nonzero()).T.tolist()
    g = ig.Graph(max(graph.shape), edgelist, directed=directed)
    # set vertices as weights
    g.es["weights"] = graph.data

    if not partition_type:
        partition_type = leidenalg.RBConfigurationVertexPartition
    if resolution_parameter:
        kargs["resolution_parameter"] = resolution_parameter
    if use_weights:
        kargs["weights"] = np.array(g.es["weights"]).astype("float64")
    kargs["n_iterations"] = n_iterations
    kargs["seed"] = seed

    print("Running Leiden optimization", flush=True)
    tic_ = time.time()
    communities = leidenalg.find_partition(
        g,
        partition_type=partition_type,
        **kargs,
    )
    Q = communities.q
    print(
        "Leiden completed in {} seconds".format(time.time() - tic_),
        flush=True,
    )
    communities = np.asarray(communities.membership)

    return communities, Q


def run_louvain(
    graph: sp.coo_matrix, q_tol: float, louvain_time_limit: int
) -> Tuple[np.ndarray, float]:
    """
    Wrapper for Louvain community detection

    Args:
        graph (sp.coo_matrix): See below in 'cluster()'
        q_tol (float): See below in 'cluster()'
        louvain_time_limit (int): See below in 'cluster()'

    Returns:
        communities, Q (Tuple[np.ndarray, float]): See below in 'cluster()'
    """
    # write to file with unique id
    uid = uuid.uuid1().hex
    graph2binary(uid, graph)
    communities, Q = runlouvain(uid, tol=q_tol, time_limit=louvain_time_limit)

    # clean up
    for f in os.listdir():
        if re.search(uid, f):
            os.remove(f)

    return communities, Q


def cluster(
    data: Union[np.ndarray, spmatrix],
    clustering_algo: Union["louvain", "leiden"] = "louvain",
    k: int = 30,
    directed: bool = False,
    prune: bool = False,
    min_cluster_size: int = 10,
    jaccard: bool = True,
    primary_metric: Union[
        "euclidean", "manhattan", "correlation", "cosine"
    ] = "euclidean",
    n_jobs: int = -1,
    q_tol: float = 1e-3,
    louvain_time_limit: int = 2000,
    nn_method: Union["kdtree", "brute"] = "kdtree",
    partition_type: Optional[Type[MutableVertexPartition]] = None,
    resolution_parameter: float = 1,
    n_iterations: int = -1,
    use_weights: bool = True,
    seed: Optional[int] = None,
    **kargs,
) -> Tuple[np.array, spmatrix, float]:
    """\
    PhenoGraph clustering

    Parameters
    ----------
    data
        Numpy ndarray of data to cluster, or sparse matrix of k-nearest neighbor graph.
        If ndarray, n-by-d array of n cells in d dimensions.
        If sparse matrix, n-by-n adjacency matrix.
    clustering_algo
        Choose `'louvain'` or `'leiden'`. Any other value will return only graph object.
    k
        Number of nearest neighbors to use in first step of graph construction.
    directed
        Whether to use a symmetric (default) or asymmetric ("directed") graph.
        The graph construction process produces a directed graph, which is symmetrized
        by one of two methods (see below).
    prune
        Whether to symmetrize by taking the average (prune = False) or product
        (prune = True) between the graph and its transpose.
    min_cluster_size
        Cells that end up in a cluster smaller than min_cluster_size are considered
        outliers and are assigned to -1 in the cluster labels.
    jaccard
        If True, use Jaccard metric between k-neighborhoods to build graph.
        If False, use a Gaussian kernel.
    primary_metric
        Distance metric to define nearest neighbors. Options include: {'euclidean',
        'manhattan', 'correlation', 'cosine'}. Note that performance will be slower for
        `correlation` and `cosine`.
    n_jobs
        Nearest Neighbors and Jaccard coefficients will be computed in parallel using
        n_jobs. If 1 is given, no parallelism is used. If set to -1, all CPUs are used.
        For n_jobs below -1, `n_cpus + 1 + n_jobs` are used.
    q_tol
        Tolerance (i.e., precision) for monitoring modularity optimization
    louvain_time_limit
        Maximum number of seconds to run modularity optimization. If exceeded the best
        result so far is returned.
    nn_method
        Whether to use brute force or kdtree for nearest neighbor search. For very large
        high-dimensional data sets, brute force (with parallel computation) performs
        faster than kdtree.
    partition_type
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`. For the
        available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    resolution_parameter
        A parameter value controlling the coarseness of the clustering in Leiden. Higher
        values lead to more clusters. Set to `None` if overriding `partition_type` to
        one that doesn’t accept a `resolution_parameter`.
    n_iterations
        Number of iterations to run the Leiden algorithm. If the number of iterations is
        negative, the Leiden algorithm is run until an iteration in which there was no
        improvement.
    use_weights
        Use vertices in the Leiden computation.
    seed
        Leiden initialization of the optimization.
    kargs
        Additional arguments passed to :func:`~leidenalg.find_partition` and the
        constructor of the `partition_type`.

    Returns
    -------
    communities
        numpy integer array of community assignments for each row in data.
    graph
        numpy sparse array of the graph that was used for clustering.
    Q
        the modularity score for communities on graph.

    Example
    -------
    >>> import phenograph
    >>> import scipy.sparse
    >>> import numpy as np

    >>> N = 5000
    >>> K = 30
    >>> RowInd = np.repeat(np.arange(N), K)
    >>> ColInd = np.tile(np.arange(N), K)
    >>> Mat = scipy.sparse.csr_matrix(
    ...     (np.ones(ColInd.shape), (RowInd, ColInd)), shape=(N, N)
    ... )

    >>> communities, graph, Q = phenograph.cluster(Mat, clustering_algo = 'leiden')
    """

    # NB if prune=True, graph must be undirected, and the prune setting takes precedence
    if prune:
        print("Setting directed=False because prune=True")
        directed = False

    if n_jobs == 1:
        kernel = jaccard_kernel
    else:
        kernel = parallel_jaccard_kernel
    kernelargs = {}

    # Start timer
    tic = time.time()
    # Go!
    if isinstance(data, sp.spmatrix) and data.shape[0] == data.shape[1]:
        print(
            "Using neighbor information from provided graph, "
            "rather than computing neighbors directly",
            flush=True,
        )
        lilmatrix = data.tolil()
        d = np.vstack(lilmatrix.data).astype("float32")  # distances
        idx = np.vstack(lilmatrix.rows).astype("int32")  # neighbor indices by row
        del lilmatrix
        assert idx.shape[0] == data.shape[0]
    else:
        d, idx = find_neighbors(
            data, k=k, metric=primary_metric, method=nn_method, n_jobs=n_jobs
        )
        print("Neighbors computed in {} seconds".format(time.time() - tic), flush=True)

    subtic = time.time()
    kernelargs["idx"] = idx
    # if not using jaccard kernel, use gaussian
    if not jaccard:
        kernelargs["d"] = d
        kernelargs["sigma"] = 1.0
        kernel = gaussian_kernel
        graph = neighbor_graph(kernel, kernelargs)
        print(
            "Gaussian kernel graph constructed in {} seconds".format(
                time.time() - subtic
            ),
            flush=True,
        )
    else:
        del d
        graph = neighbor_graph(kernel, kernelargs)
        print(
            "Jaccard graph constructed in {} seconds".format(time.time() - subtic),
            flush=True,
        )
    if not directed:
        if not prune:
            # symmetrize graph by averaging with transpose
            sg = (graph + graph.transpose()).multiply(0.5)
        else:
            # symmetrize graph by multiplying with transpose
            sg = graph.multiply(graph.transpose())
        # retain lower triangle (for efficiency)
        graph = sp.tril(sg, -1)

    # choose between Louvain or Leiden algorithm
    communities, Q = "", ""
    if clustering_algo == "louvain":
        communities, Q = run_louvain(graph, q_tol, louvain_time_limit)

    elif clustering_algo == "leiden":
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

    else:
        # return only graph object
        pass

    print("Sorting communities by size, please wait ...", flush=True)
    communities = sort_by_size(communities, min_cluster_size)

    print("PhenoGraph completed in {} seconds".format(time.time() - tic), flush=True)

    return communities, graph, Q
