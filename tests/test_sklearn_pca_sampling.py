import numpy as np

from genomics.core.sklearn_pca_cache import _stratified_fit_indices


def test_stratified_fit_indices_preserve_classes_and_size():
    labels = np.asarray([0] * 100 + [1] * 50 + [2] * 10)

    indices = _stratified_fit_indices(
        labels,
        fraction=0.5,
        min_samples=90,
        min_samples_per_class=5,
        random_seed=13,
        stratify=True,
    )

    selected_labels = labels[indices]
    assert len(indices) == 90
    assert set(selected_labels.tolist()) == {0, 1, 2}
    assert int((selected_labels == 2).sum()) >= 5


def test_stratified_fit_indices_returns_all_when_request_exceeds_dataset():
    labels = np.asarray([0, 0, 1, 1])

    indices = _stratified_fit_indices(
        labels,
        fraction=0.5,
        min_samples=10,
        min_samples_per_class=1,
        random_seed=13,
        stratify=True,
    )

    assert indices == [0, 1, 2, 3]


def test_stratified_fit_indices_can_use_train_only_subset():
    labels = np.asarray([0] * 40 + [1] * 40)

    indices = _stratified_fit_indices(
        labels,
        fraction=0.25,
        min_samples=30,
        min_samples_per_class=10,
        random_seed=7,
        stratify=True,
    )

    selected_labels = labels[indices]
    assert len(indices) == 30
    assert int((selected_labels == 0).sum()) >= 10
    assert int((selected_labels == 1).sum()) >= 10
