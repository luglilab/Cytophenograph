from phenograph.classify import classify, random_walk_probabilities


def test_classify_generated(cluster_dataset_generated):
    output, _ = classify(
        cluster_dataset_generated.train_set, cluster_dataset_generated.test_set
    )
    actual = cluster_dataset_generated.test_labels
    assert all([p == a for p, a in zip(output, actual)])


def test_classify_fixed(cluster_dataset_fixed):
    output, _ = classify(
        cluster_dataset_fixed.train_set, cluster_dataset_fixed.test_set, k=5
    )
    actual = cluster_dataset_fixed.actual_labels
    assert all([p == a for p, a in zip(output, actual)])


def test_random_walk_probabilities_fixed(cluster_dataset_fixed):
    from numpy.testing import assert_array_equal

    output = random_walk_probabilities(
        cluster_dataset_fixed.sparse_affinity_matrix,
        cluster_dataset_fixed.initial_labels,
    )
    actual = cluster_dataset_fixed.random_walk_P
    assert_array_equal(output, actual)
