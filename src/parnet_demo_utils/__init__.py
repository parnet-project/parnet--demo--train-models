"""parnet_demo_utils — local utilities for PARNET eCLIP demo notebooks.

Provides dataset filtering, sparse tensor helpers, HFDS conversion utilities,
and BED file parsing for the notebooks in this repository.

Quick start::

    from functools import partial
    from parnet_demo_utils import (
        FilteredMultiTaskDataset,
        filter_min_read_count,
        center_crop_item,
        parse_tile_name,
    )

    filters = [
        partial(filter_min_read_count,
                min_read_count=3,
                tasks=["eCLIP"],
                track_indices=my_indices),
    ]
    ds = FilteredMultiTaskDataset(data["train"], my_indices, filters)
    sample = ds[0]
"""

from .bed_utils import GenomicInterval, parse_tile_name
from .datasets import FilteredMultiTaskDataset, center_crop_item
from .extracted_datasets import prepare_chr8_gff, prepare_chr8_peaks, prepare_chr8_splice_sites
from .filters import FilterFunction, filter_min_read_count, filter_minimum_length
from .hfds_utils import convert_hfds_sample, hfds_sample_passes_filter
from .splice_sites import extract_splice_sites_from_gff, extract_splice_sites_from_transcript_gff
from .training_utils import MetricHistory
from .sparse_utils import (
    ParnetDataElement,
    ParnetDataElementInputs,
    ParnetDataElementOutputs,
    SparseTensorDict,
    classify_seq_padding,
    dense_onehot_to_seq,
    dense_to_sparse_lists,
    get_padding_mask,
    infer_pad_sizes,
    sparse_onehot_to_seq,
    strip_element_padding,
    torch_dense_to_sparse,
    torch_sparse_to_dense,
)

__all__ = [
    # sparse_utils
    "SparseTensorDict",
    "ParnetDataElementInputs",
    "ParnetDataElementOutputs",
    "ParnetDataElement",
    "torch_sparse_to_dense",
    "torch_dense_to_sparse",
    "dense_to_sparse_lists",
    "sparse_onehot_to_seq",
    "dense_onehot_to_seq",
    "classify_seq_padding",
    "infer_pad_sizes",
    "get_padding_mask",
    "strip_element_padding",
    # filters
    "FilterFunction",
    "filter_minimum_length",
    "filter_min_read_count",
    # datasets
    "FilteredMultiTaskDataset",
    "center_crop_item",
    # hfds_utils
    "hfds_sample_passes_filter",
    "convert_hfds_sample",
    # bed_utils
    "GenomicInterval",
    "parse_tile_name",
    # splice_sites
    "extract_splice_sites_from_transcript_gff",
    "extract_splice_sites_from_gff",
    # extracted_datasets
    "prepare_chr8_gff",
    "prepare_chr8_splice_sites",
    "prepare_chr8_peaks",
    # training_utils
    "MetricHistory",
]
