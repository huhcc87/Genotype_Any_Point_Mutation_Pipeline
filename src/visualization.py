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
    """Schematic gel with actual fragment-size labels beside each predicted lane."""
    if ladder_key not in LADDERS:
        raise ValueError("Unknown ladder preset.")
    marker = LADDERS[ladder_key]
    profiles = gel_fragment_profiles(amplicon_size, enzyme["cuts_wt"], enzyme["cuts_mut"])
    labels = ["DNA ladder"] + list(profiles)
    lanes = [marker["sizes"]] + list(profiles.values())
    sizes = [int(s) for lane in lanes for s in lane if s > 0]
    if not sizes:
        raise ValueError("No predicted fragments.")
    upper = max(400, int(max(sizes) * 1.32))
    lower = 18
    def y_pos(bp):
        return .14 + .72 * ((math.log(upper) - math.log(max(bp, lower))) /
                             (math.log(upper) - math.log(lower)))

    fig, ax = plt.subplots(figsize=(14.5, 8.8))
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#142238")
    ax.set_xlim(.18, 6.15)
    ax.set_ylim(1.04, .025)
    ax.set_xticks(range(1, 6), labels, rotation=10, ha="right",
                  color="#1b293c", fontsize=10)
    ax.tick_params(axis="x", length=0, pad=10)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Distinguish fragments close enough to overlap by placing label strings
    # in a right-side summary table below the graph, not by altering sizes.
    for lane_index, lane in enumerate(lanes, 1):
        ax.add_patch(Rectangle((lane_index-.23, .092), .46, .018,
                               color="#a2b3c8"))
        for size,count in Counter(int(x) for x in lane if x>=lower).items():
            y=y_pos(size)
            ax.add_patch(Rectangle((lane_index-.24,y-.0048),.48,.0096,
                                   color=(.24,.94,.66,min(.98,.55+.1*count))))
            if lane_index == 1:
                ax.text(.64,y,f"{size:,} bp",color="#f0f7ff",fontsize=8,
                        va="center",ha="right")
            else:
                ax.text(lane_index+.27,y,f"{size:,}",color="#f0f7ff",
                        fontsize=8,va="center",ha="left")
    ax.set_title(f"Predicted complete digest: {enzyme['enzyme']}",
                 color="#ffffff",fontsize=14,pad=14)
    fig.text(.12,.14,"Sample fragments (bp):  "+"   |   ".join(
        f"{k}: {', '.join(map(str,sorted(v,reverse=True)))}"
        for k,v in profiles.items()),fontsize=9,color="#1e293b",wrap=True)
    fig.text(.12,.105,"Selected marker: "+ladder_key,
             fontsize=9,color="#334155")
    fig.text(.12,.073,ladder_notes(ladder_key,[x for v in profiles.values() for x in v]),
             fontsize=8,color="#475569",wrap=True)
    fig.text(.12,.04,"Simulation only; shared-size DNA bands co-migrate. "
             "Real gel migration and intensity are not calibrated.",fontsize=8,color="#64748b")
    fig.subplots_adjust(left=.08,right=.96,top=.89,bottom=.23)
    target=Path(output_path)
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target,dpi=190,facecolor=fig.get_facecolor())
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
