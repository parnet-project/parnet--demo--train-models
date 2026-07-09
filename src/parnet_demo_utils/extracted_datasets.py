"""Idempotent preparation of chr8-subsetted demo resources.

Both ``embeddings.py.ipynb`` and ``intepretation.py.ipynb`` under
``notebooks/demo--train-spliceosome-hepg2/`` need a small, fast-to-load slice
of two large NAS-hosted resources:

- The MANE GFF (genome-wide annotation), subset to chromosome 8 — one of the
  two held-out test chromosomes for this model — from which donor/acceptor
  splice sites are extracted.
- PRPF8/U2AF2 HepG2 eCLIP IDR narrowPeak files, subset to chromosome 8.

Each ``prepare_*`` function writes into
``results/spliceosome-hepg2.precomputed/extracted_datasets/`` and skips work
if its output already exists, so either notebook can call them independently
and re-running is cheap.
"""

from __future__ import annotations

import datetime
import gzip
import os
from pathlib import Path

import pandas as pd
import yaml
from pylbsr.bio.gff import ExtendedGFF

from .splice_sites import extract_splice_sites_from_gff

NARROWPEAK_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
]


def _update_metadata(output_dir: os.PathLike, key: str, entry: dict) -> None:
    """Merge ``entry`` under ``key`` into ``output_dir/metadata.yaml``."""
    output_dir = Path(output_dir)
    metadata_path = output_dir / "metadata.yaml"
    metadata = yaml.safe_load(metadata_path.read_text()) if metadata_path.exists() else {}
    metadata[key] = {**entry, "generated_at": datetime.date.today().isoformat()}
    metadata_path.write_text(yaml.dump(metadata, sort_keys=False))


def prepare_chr8_gff(mane_gff_path: os.PathLike, output_dir: os.PathLike, chrom: str = "chr8") -> Path:
    """Subset a MANE GFF file to a single chromosome by a plain line filter.

    Filters by the ``seqid`` column (first tab-delimited field) directly on
    the raw text, without loading the (genome-wide, ~525k-line) file into
    pandas — only the much smaller chromosome subset is parsed downstream.

    Args:
        mane_gff_path: Path to the full MANE ``*.ensembl_genomic.gff`` file.
        output_dir: Directory to write ``{chrom}.MANE.ensembl_genomic.gff``
            into (created if missing).
        chrom: Chromosome to keep, e.g. ``"chr8"``.

    Returns:
        Path to the written (or pre-existing) subsetted GFF file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{chrom}.MANE.ensembl_genomic.gff"

    if output_path.exists():
        return output_path

    n_rows = 0
    with open(mane_gff_path) as f_in, open(output_path, "w") as f_out:
        f_out.write("##gff-version 3\n")
        for line in f_in:
            if line.startswith("#"):
                continue
            if line.split("\t", 1)[0] == chrom:
                f_out.write(line)
                n_rows += 1

    _update_metadata(
        output_dir,
        "chr8_gff",
        {"source": str(mane_gff_path), "chrom": chrom, "n_rows": n_rows, "output": output_path.name},
    )
    return output_path


def prepare_chr8_splice_sites(chr8_gff_path: os.PathLike, output_dir: os.PathLike) -> Path:
    """Extract donor/acceptor splice sites from a chromosome-subsetted GFF.

    Args:
        chr8_gff_path: Path to a chr8-subsetted GFF (as produced by
            :func:`prepare_chr8_gff`).
        output_dir: Directory to write
            ``chr8.splice_sites.donor_acceptor.tsv`` into.

    Returns:
        Path to the written (or pre-existing) splice-sites TSV, with columns
        ``chrom, ss_pos, strand, ss_type, transcript_id, gene_name``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "chr8.splice_sites.donor_acceptor.tsv"

    if output_path.exists():
        return output_path

    gff_extended = ExtendedGFF.from_filepath(chr8_gff_path).extended
    splice_sites_df = extract_splice_sites_from_gff(gff_extended)
    splice_sites_df.to_csv(output_path, sep="\t", index=False)

    _update_metadata(
        output_dir,
        "chr8_splice_sites",
        {
            "source": str(chr8_gff_path),
            "n_sites": len(splice_sites_df),
            "n_donor": int((splice_sites_df["ss_type"] == "donor").sum()),
            "n_acceptor": int((splice_sites_df["ss_type"] == "acceptor").sum()),
            "output": output_path.name,
        },
    )
    return output_path


def prepare_chr8_peaks(rbp_ct: str, narrowpeak_path: os.PathLike, output_dir: os.PathLike, chrom: str = "chr8") -> Path:
    """Subset a narrowPeak eCLIP peak file to a single chromosome.

    Args:
        rbp_ct: RBP + cell-type identifier used in the output filename, e.g.
            ``"PRPF8_HepG2"``.
        narrowpeak_path: Path to the source ``narrowPeak.*.bed.gz`` file
            (bed10 format: chrom, start, end, name, score, strand,
            signalValue, pValue, qValue, peak).
        output_dir: Directory to write ``{rbp_ct}.chr8.narrowPeak.bed.gz``
            into.
        chrom: Chromosome to keep, e.g. ``"chr8"``.

    Returns:
        Path to the written (or pre-existing) subsetted, gzipped BED file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{rbp_ct}.chr8.narrowPeak.bed.gz"

    if output_path.exists():
        return output_path

    peaks = pd.read_csv(narrowpeak_path, sep="\t", header=None, names=NARROWPEAK_COLUMNS)
    peaks_chrom = peaks.loc[peaks["chrom"] == chrom].reset_index(drop=True)

    with gzip.open(output_path, "wt") as f_out:
        peaks_chrom.to_csv(f_out, sep="\t", header=False, index=False)

    _update_metadata(
        output_dir,
        f"{rbp_ct}_chr8_peaks",
        {"source": str(narrowpeak_path), "chrom": chrom, "n_peaks": len(peaks_chrom), "output": output_path.name},
    )
    return output_path
