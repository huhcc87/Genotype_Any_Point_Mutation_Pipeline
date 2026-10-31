"""Transcript-specific gene architecture from verified reference exon coordinates.

Coordinates displayed to users are 1-based inclusive. Exons include UTRs,
unlike the app's separate coding-exon view. Intron counts apply to ONE chosen
transcript, not all alternative transcript isoforms of a gene.
"""
import json


def _length(start, end):
    return end - start + 1


def gene_locus_summary(gene):
    """Gene-level locus metadata is optional and never inferred from one isoform."""
    chromosome = gene.get("chromosome")
    locus = gene.get("gene_locus") or {}
    start = locus.get("start_1based")
    end = locus.get("end_1based")
    if start is not None or end is not None:
        if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
            raise ValueError("Gene locus coordinates are invalid.")
    return {
        "gene_symbol": gene.get("gene_symbol"),
        "gene_id": gene.get("gene_id"),
        "assembly": gene.get("assembly"),
        "chromosome": chromosome,
        "cytogenetic_band": gene.get("cytoband") or None,
        "gene_start_1based": start,
        "gene_end_1based": end,
        "gene_span_bp": _length(start, end) if start is not None else None,
        "gene_span_kb": round(_length(start, end) / 1000, 3) if start is not None else None,
        "gene_locus_source": gene.get("gene_locus_source") or "Not available",
        "scope": "Gene locus (all isoforms) if provided by the source; not selected-transcript coordinates.",
    }


def _outside_cds(start, end, cds_low, cds_high, strand):
    before = max(0, min(end, cds_low - 1) - start + 1)
    after = max(0, end - max(start, cds_high + 1) + 1)
    coding = max(0, min(end, cds_high) - max(start, cds_low) + 1)
    # On the negative strand, lower genomic coordinates are 3' rather than 5'.
    utr5, utr3 = (before, after) if strand == 1 else (after, before)
    return coding, utr5, utr3


def transcript_structure(transcript, cds=None, codon_rows=None):
    """Return summary + exon/intron/junction tables for a selected transcript.

    Exons are sorted in TRANSCRIPT direction, with genomic start/end always
    reported in ascending genomic order. Junction donor/acceptor positions
    are theoretical intronic windows; canonical GT/AG motifs are not asserted.
    """
    exons = transcript.get("Exon") or []
    if not isinstance(exons, list) or not exons:
        raise ValueError("No exon model available for selected transcript.")
    strand = transcript.get("strand")
    if strand not in (1, -1):
        raise ValueError("Transcript strand must be +1 or -1.")
    normalized = []
    for e in exons:
        low, high = int(e["start"]), int(e["end"])
        if low < 1 or high < low:
            raise ValueError("Invalid exon coordinates.")
        normalized.append((low, high))
    ordered = sorted(normalized, reverse=strand == -1)
    genomic_order = sorted(normalized)
    if any(genomic_order[i][1] >= genomic_order[i+1][0] for i in range(len(genomic_order)-1)):
        raise ValueError("Exon intervals overlap; cannot calculate introns reliably.")

    tx_start, tx_end = genomic_order[0][0], genomic_order[-1][1]
    translation = transcript.get("Translation") or {}
    if translation:
        coding_low = min(int(translation["start"]),int(translation["end"]))
        coding_high = max(int(translation["start"]),int(translation["end"]))
        if coding_low < tx_start or coding_high > tx_end:
            raise ValueError("CDS boundaries are outside the transcript.")
    else:
        coding_low = coding_high = None

    exon_table = []
    spliced_offset = 1
    for number, (low, high) in enumerate(ordered, 1):
        coding, utr5, utr3 = (
            _outside_cds(low, high, coding_low, coding_high, strand)
            if coding_low is not None else (0, None, None)
        )
        length = _length(low, high)
        if coding_low is not None and coding + utr5 + utr3 != length:
            raise ValueError("Transcript exon/CDS lengths do not reconcile.")
        region = ("coding" if coding == length else "UTR-only" if coding == 0
                  else "mixed CDS/UTR") if coding_low is not None else "CDS not available"
        exon_table.append({
            "exon_number": number,
            "genomic_start_1based": low,
            "genomic_end_1based": high,
            "exon_length_bp": length,
            "transcript_start_1based": spliced_offset,
            "transcript_end_1based": spliced_offset + length - 1,
            "coding_bases": coding,
            "utr_5prime_bp": utr5,
            "utr_3prime_bp": utr3,
            "region_type": region,
            "genomic_5prime_end": high if strand == -1 else low,
            "genomic_3prime_end": low if strand == -1 else high,
        })
        spliced_offset += length

    intron_table = []
    junction_table = []
    for i, (upstream, downstream) in enumerate(zip(ordered[:-1], ordered[1:]), 1):
        if strand == 1:
            lower, upper = upstream[1] + 1, downstream[0] - 1
            donor = [upstream[1] + 1, upstream[1] + 2]
            acceptor = [downstream[0] - 2, downstream[0] - 1]
            donor_exon_edge, acceptor_exon_edge = upstream[1], downstream[0]
        else:
            lower, upper = downstream[1] + 1, upstream[0] - 1
            donor = [upstream[0] - 1, upstream[0] - 2]
            acceptor = [downstream[1] + 2, downstream[1] + 1]
            donor_exon_edge, acceptor_exon_edge = upstream[0], downstream[1]
        if lower > upper:
            raise ValueError("Adjacent exons lack a valid intervening intron.")
        length = _length(lower, upper)
        intron_table.append({
            "intron_number": i, "between_exons": f"{i}–{i+1}",
            "genomic_start_1based":lower,"genomic_end_1based":upper,
            "intron_length_bp":length,"strand":"+" if strand == 1 else "-",
        })
        junction_table.append({
            "junction_number":i,
            "upstream_exon":i,"downstream_exon":i+1,
            "donor_exon_boundary_1based":donor_exon_edge,
            "acceptor_exon_boundary_1based":acceptor_exon_edge,
            "donor_intronic_positions_5to3":donor if length >= 4 else None,
            "acceptor_intronic_positions_5to3":acceptor if length >= 4 else None,
            "intron_length_bp":length,
            "motifs_verified":False,
        })
    spliced_length = sum(e["exon_length_bp"] for e in exon_table)
    total_intronic = sum(i["intron_length_bp"] for i in intron_table)
    if spliced_length + total_intronic != _length(tx_start, tx_end):
        raise ValueError("Exons and introns do not add up to transcript genomic span.")

    cds_length = len(cds) if cds is not None else (
        sum(e["coding_bases"] for e in exon_table) if coding_low is not None else None
    )
    if cds_length is not None and coding_low is not None and cds_length != sum(e["coding_bases"] for e in exon_table):
        raise ValueError("CDS sequence length and annotated coding exon length disagree.")
    terminal_stop = bool(cds and cds[-3:] in ("TAA","TAG","TGA"))
    protein_length = (len(cds)//3 - int(terminal_stop)) if cds is not None else None
    summary = {
        "transcript_id": transcript.get("id"),
        "chromosome": transcript.get("seq_region_name"),
        "strand": "+" if strand == 1 else "-",
        "transcript_genomic_start_1based":tx_start,
        "transcript_genomic_end_1based":tx_end,
        "transcript_genomic_span_bp":_length(tx_start,tx_end),
        "spliced_transcript_length_bp":spliced_length,
        "total_exons":len(exon_table),
        "total_introns":len(intron_table),
        "coding_exons":sum(e["coding_bases"]>0 for e in exon_table),
        "noncoding_exons":sum(e["coding_bases"]==0 for e in exon_table) if coding_low is not None else None,
        "cds_length_bp":cds_length,
        "coding_triplets_including_stop":len(cds)//3 if cds is not None else None,
        "terminal_stop_codon":cds[-3:] if terminal_stop else None,
        "predicted_protein_length_aa":protein_length,
        "utr_5prime_total_bp":sum(e["utr_5prime_bp"] for e in exon_table) if coding_low is not None else None,
        "utr_3prime_total_bp":sum(e["utr_3prime_bp"] for e in exon_table) if coding_low is not None else None,
        "intronic_total_bp":total_intronic,
        "coordinate_convention":"1-based inclusive genomic coordinates; exon numbering in transcript 5'→3' order",
        "reference_scope":"Selected transcript model only; other isoforms can differ.",
    }
    if codon_rows is not None:
        if cds is None or len(codon_rows) != len(cds)//3:
            raise ValueError("Codon table and selected CDS do not reconcile.")
        stop_rows = [r for r in codon_rows if r["reference_amino_acid"] == "*"]
        summary["stop_codon_positions"] = [r["codon_number"] for r in stop_rows]
        summary["internal_stop_warning"] = any(r["codon_number"] < len(codon_rows) for r in stop_rows)
    return {"summary":summary,"exons":exon_table,"introns":intron_table,"junctions":junction_table}
