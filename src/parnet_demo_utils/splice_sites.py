"""Strand-aware donor/acceptor splice site extraction from GFF annotations.

Ports the extraction logic validated in
``parnet--analyses--mutations/src/mutations_analyses_libs/splice_sites.py``
so this repo does not need to depend on that sibling project — everything
here is built on :mod:`pylbsr.bio.gff`, which is already a dependency.
"""

from __future__ import annotations

import pandas as pd

from pylbsr.bio.gff import gff_transcript_segments_to_bed


def extract_splice_sites_from_transcript_gff(gff_transcript: pd.DataFrame) -> list[dict]:
    """Identify donor/acceptor splice site positions for a single transcript.

    Uses :func:`pylbsr.bio.gff.gff_transcript_segments_to_bed` to segment the
    transcript into exon + intron BED intervals (always sorted by ascending
    genomic coordinate), then labels each exon-intron boundary as a donor or
    acceptor site. No GT-AG dinucleotide restriction is applied.

    Because segments are sorted by ascending genomic coordinate, the first
    intronic base (genomically) is a donor on the ``+`` strand but an
    acceptor on the ``-`` strand (transcription runs right-to-left there);
    the last intronic base is the reverse.

    Args:
        gff_transcript: GFF dataframe for exactly one transcript (must have
            ``transcript_id``, ``ID``, ``type`` columns and a single
            ``transcript`` row plus at least one ``exon`` row). Rows of other
            types (CDS, UTR, ...) are dropped before segmentation so they
            don't interfere with intron detection.

    Returns:
        List of dicts with keys ``chrom``, ``ss_pos`` (1-based genomic
        position), ``strand``, ``ss_type`` (``"donor"`` or ``"acceptor"``),
        ``transcript_id``. Empty if the transcript has no exon or no
        transcript row.
    """
    gff_for_seg = gff_transcript.loc[gff_transcript["type"].isin(["transcript", "exon"])].copy()

    if "transcript" not in gff_for_seg["type"].values:
        return []
    if "exon" not in gff_for_seg["type"].values:
        return []

    bed = gff_transcript_segments_to_bed(gff_for_seg)

    strand = bed["strand"].iloc[0]
    chrom = bed["chrom"].iloc[0]
    transcript_id = gff_transcript["transcript_id"].iloc[0]

    splice_sites = []
    for i in range(len(bed) - 1):
        seg_a = bed.iloc[i]
        seg_b = bed.iloc[i + 1]

        a_is_intron = str(seg_a["name"]).startswith("intron:")
        b_is_intron = str(seg_b["name"]).startswith("intron:")

        # Exon -> intron boundary: first intronic base (genomically).
        if (not a_is_intron) and b_is_intron:
            ss_pos = int(seg_b["start"]) + 1  # BED 0-based start -> 1-based
            splice_sites.append(
                dict(
                    chrom=chrom,
                    ss_pos=ss_pos,
                    strand=strand,
                    ss_type="donor" if strand != "-" else "acceptor",
                    transcript_id=transcript_id,
                )
            )

        # Intron -> exon boundary: last intronic base (genomically).
        elif a_is_intron and (not b_is_intron):
            ss_pos = int(seg_a["end"])  # BED 0-based exclusive end -> 1-based inclusive
            splice_sites.append(
                dict(
                    chrom=chrom,
                    ss_pos=ss_pos,
                    strand=strand,
                    ss_type="acceptor" if strand != "-" else "donor",
                    transcript_id=transcript_id,
                )
            )

    return splice_sites


def extract_splice_sites_from_gff(gff_extended: pd.DataFrame) -> pd.DataFrame:
    """Extract donor/acceptor splice sites for every transcript in a GFF.

    Args:
        gff_extended: GFF dataframe with split attribute columns (as produced
            by ``pylbsr.bio.gff.ExtendedGFF.from_filepath(...).extended``),
            covering any number of transcripts.

    Returns:
        DataFrame with columns ``chrom``, ``ss_pos``, ``strand``, ``ss_type``,
        ``transcript_id``, ``gene_name``, one row per donor/acceptor site.
        Deduplicated by ``(chrom, ss_pos, strand, ss_type)`` — shared splice
        sites across isoforms of the same gene are kept once.
    """
    gene_name_by_transcript = (
        gff_extended.loc[gff_extended["type"] == "transcript"]
        .set_index("transcript_id")["gene_name"]
        .to_dict()
    )

    all_splice_sites: list[dict] = []
    for transcript_id, transcript_gff in gff_extended.groupby("transcript_id", sort=False):
        if pd.isna(transcript_id):
            continue
        splice_sites = extract_splice_sites_from_transcript_gff(transcript_gff)
        for ss in splice_sites:
            ss["gene_name"] = gene_name_by_transcript.get(transcript_id)
        all_splice_sites.extend(splice_sites)

    splice_sites_df = pd.DataFrame(
        all_splice_sites,
        columns=["chrom", "ss_pos", "strand", "ss_type", "transcript_id", "gene_name"],
    )
    return splice_sites_df.drop_duplicates(
        subset=["chrom", "ss_pos", "strand", "ss_type"]
    ).reset_index(drop=True)
