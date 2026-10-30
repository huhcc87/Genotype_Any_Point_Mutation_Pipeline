"""Primer-position annotations and schematic, size-labelled PCR-RFLP gels.

Figures are in-silico *illustrations*. They do not calibrate gel migration,
fragment brightness, digestion efficiency, or marker mass.
"""
from collections import Counter
from html import escape
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from .gel import gel_fragment_profiles
from .ladders import LADDERS, ladder_notes


def render_ladder_gel(amplicon_size, enzyme, ladder_key, output_path):
    """Plot the DNA ladder and WT/heterozygous/mutant digest lanes."""
    if ladder_key not in LADDERS:
        raise ValueError("Unknown ladder preset.")
    marker = LADDERS[ladder_key]
    profiles = gel_fragment_profiles(amplicon_size, enzyme["cuts_wt"], enzyme["cuts_mut"])
    labels = ["DNA ladder"] + list(profiles)
    lanes = [marker["sizes"]] + list(profiles.values())
    all_sizes = [s for group in lanes for s in group if s > 0]
    if not all_sizes:
        raise ValueError("No DNA fragments to display.")
    upper = max(500, int(max(all_sizes) * 1.25))
    lower = 18

    def y_position(size):
        return 0.1 + 0.78 * (
            (math.log(upper) - math.log(max(lower, size))) /
            (math.log(upper) - math.log(lower))
        )

    fig, ax = plt.subplots(figsize=(11, 6.4))
    fig.patch.set_facecolor("#f9fbfd")
    ax.set_facecolor("#152238")
    ax.set_xlim(0.45, len(labels) + 0.55)
    ax.set_ylim(1.02, 0.02)

    for i, fragments in enumerate(lanes, 1):
        ax.add_patch(Rectangle((i - .27, .062), .54, .018, color="#96a9bd"))
        counts = Counter(s for s in fragments if s >= lower)
        for size, repeats in counts.items():
            ax.add_patch(Rectangle(
                (i - .27, y_position(size) - .006), .54, .012,
                color=(.25, .92, .66, min(.95, .49 + .13 * repeats)),
            ))

    ax.set_xticks(range(1, len(labels) + 1), labels, rotation=10, ha="right")
    ax.tick_params(axis="x", colors="#27394c", length=0, pad=10)
    ticks = [t for t in [25, 50, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000,
                         3000, 5000, 10000] if lower <= t <= upper]
    ax.set_yticks([y_position(t) for t in ticks], [f"{t:,}" for t in ticks])
    ax.tick_params(axis="y", colors="white", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        f"Predicted complete digest: {enzyme['enzyme']} (fragment size, bp)",
        color="#ffffff", pad=14
    )
    description = " | ".join(
        f"{label}: {', '.join(map(str, sorted(fragments, reverse=True)))} bp"
        for label, fragments in profiles.items()
    )
    fig.text(.10, .078, description, fontsize=8, color="#334155", wrap=True)
    footnote = f"Ladder: {ladder_key}. " + ladder_notes(
        ladder_key, [x for group in profiles.values() for x in group]
    )
    fig.text(.10, .041, footnote, fontsize=7.5, color="#475569", wrap=True)
    fig.subplots_adjust(left=.10, right=.98, top=.87, bottom=.29)

    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    return str(target)


def sequence_markup(result, max_bases=1600):
    """Self-contained color-highlighted plus-strand amplicon HTML."""
    pair = result.get("selected_primer")
    dna = result.get("wt_amplicon")
    if not pair or not dna:
        return "<p>No designed amplicon is available.</p>"

    start = pair["amplicon_start_1based"]
    variant_index = result["position"] - start
    forward_end = len(pair["forward"])
    reverse_start = len(dna) - len(pair["reverse"])
    left = max(0, variant_index - max_bases // 2)
    right = min(len(dna), left + max_bases)
    left = max(0, right - max_bases)
    rows = []
    for offset in range(left, right, 60):
        characters = []
        for i in range(offset, min(offset + 60, right)):
            css = ("mutation" if i == variant_index else
                   "forward" if i < forward_end else
                   "reverse" if i >= reverse_start else "other")
            characters.append(f'<span class="{css}">{escape(dna[i])}</span>')
        rows.append(
            f'<div><span class="position">{start + offset:>9,}</span> '
            + "".join(characters) + "</div>"
        )

    html = """<!DOCTYPE html><html><head><meta charset="utf-8">
    <style>body{font:14px system-ui,Arial,sans-serif;color:#1e293b;margin:12px}
    pre{font:14px/1.8 ui-monospace,monospace;overflow-x:auto;background:#f6f9fc;
    padding:12px;border:1px solid #dae5ee;border-radius:8px}
    .position{color:#64748b;margin-right:10px}
    .forward{background:#bbf7d0;color:#14532d}
    .reverse{background:#bfdbfe;color:#1e40af}
    .mutation{background:#fde68a;color:#7c2d12;font-weight:800}
    .tag{display:inline-block;padding:4px 7px;margin-right:8px}
    </style></head><body><h3>Reference DNA and primer binding sites</h3>"""
    html += (
        f"<p>Amplicon: input positions {start:,}–{pair['amplicon_end_1based']:,} "
        f"({len(dna):,} bp); variant {escape(result['mutation'])} at "
        f"position {result['position']:,}.</p>"
    )
    html += (
        '<p><span class="tag forward">Forward binding site</span>'
        '<span class="tag mutation">WT variant base</span>'
        '<span class="tag reverse">Reverse binding site on plus strand</span></p>'
    )
    html += "<pre>" + "\n".join(rows) + "</pre>"
    html += (
        f"<p><b>Forward 5′→3′:</b> {escape(pair['forward'])} | "
        f"<b>Reverse 5′→3′:</b> {escape(pair['reverse'])}</p>"
    )
    if left or right < len(dna):
        html += "<p>Only the region around the variant is shown; full sequence is exported.</p>"
    return html + "<p>Diagram is schematic; primer specificity requires external testing.</p></body></html>"


def render_sequence_figure(result, output_path):
    """Draw primer-binding intervals and variant position on local coordinates."""
    pair = result.get("selected_primer")
    if not pair:
        raise ValueError("No selected primer pair is available.")
    start = pair["amplicon_start_1based"]
    end = pair["amplicon_end_1based"]
    mutation = result["position"]
    fwd_end = start + len(pair["forward"]) - 1
    rev_start = end - len(pair["reverse"]) + 1
    fig, ax = plt.subplots(figsize=(10, 2.6))
    fig.patch.set_facecolor("#ffffff")
    ax.plot([start, end], [0, 0], lw=5, color="#475569")
    ax.plot([start, fwd_end], [0, 0], lw=9, color="#16a34a")
    ax.plot([rev_start, end], [0, 0], lw=9, color="#2563eb")
    ax.scatter([mutation], [0], s=120, marker="D", color="#eab308", zorder=5)
    ax.text(mutation, .16, f"{result['mutation']}\n{mutation:,} bp",
            ha="center", fontsize=9)
    ax.text(start, -.23, "Forward 5′→3′\n" + pair["forward"],
            ha="left", fontsize=8, color="#166534")
    ax.text(end, -.23, "Reverse 5′→3′\n" + pair["reverse"],
            ha="right", fontsize=8, color="#1d4ed8")
    ax.set_xlim(start - 5, end + 5)
    ax.set_ylim(-.43, .41)
    ax.axis("off")
    ax.set_title(f"Primer map — {pair['product_size']} bp amplicon")
    fig.tight_layout()
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(target)
