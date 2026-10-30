"""Publication-oriented assay diagrams and manufacturer-specific ladder visuals."""
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
    """Schematic complete-digestion model, not measured electrophoretic migration."""
    marker = LADDERS[ladder_key]
    profiles = gel_fragment_profiles(
        amplicon_size, enzyme["cuts_wt"], enzyme["cuts_mut"]
    )
    names = ["DNA ladder", "WT undigested", "WT digested", "Heterozygote", "Mutant digested"]
    bands = [marker["sizes"]] + list(profiles.values())
    high = max([amplicon_size * 1.25] + [x for lane in bands for x in lane])
    low = min(20, min(x for lane in bands for x in lane))
    high = max(500, high)
    low = max(5, low)

    def ordinate(bp):
        return 0.1 + 0.8 * ((math.log10(high) - math.log10(bp)) /
                            (math.log10(high) - math.log10(low)))

    fig, ax = plt.subplots(figsize=(11.5, 7))
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#152238")
    ax.set_xlim(0.5, 5.6)
    ax.set_ylim(1.02, 0.02)
    for x in range(1, 6):
        ax.add_patch(Rectangle((x - .31, .066), .62, .018, color="#9aaac3"))
    for x, size_list in enumerate(bands, 1):
        counts = Counter(s for s in size_list if s >= 20)
        for size, copies in counts.items():
            y = ordinate(size)
            mass = (size / max(amplicon_size, 1))
            alpha = min(1., .48 + .15 * copies + .18 * mass)
            ax.add_patch(Rectangle((x-.31, y-.0055), .62, .011,
                                   color=(.25, .94, .70, alpha)))
    ax.set_xticks(range(1, 6), names, color="#263346", rotation=12, ha="right")
    ticks = sorted({t for t in ([25,50,75,100,150,200,300,400,500,750,1000,
                                 1500,2000,3000,5000,10000]) if low <= t <= high}, reverse=True)
    ax.set_yticks([ordinate(t) for t in ticks])
    ax.set_yticklabels([f"{t:,}" for t in ticks], color="white")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0, pad=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f"In-silico complete digest — {enzyme['enzyme']}  |  sizes in bp",
                 color="#e5edf8", fontsize=13, pad=18)
    fig.text(.13, .08, "Expected bands (bp): " +
             "    |    ".join(f"{k}: {','.join(map(str, sorted(v, reverse=True)))}"
                            for k, v in profiles.items()),
             ha="left", fontsize=8, color="#334155", wrap=True)
    fig.text(.13, .045, f"Ladder: {ladder_key}. " +
             ladder_notes(ladder_key, [x for lane in profiles.values() for x in lane]),
             fontsize=8, color="#334155", wrap=True)
    fig.subplots_adjust(left=.12, right=.96, top=.88, bottom=.25)
    path=Path(output_path)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=170,facecolor=fig.get_facecolor())
    plt.close(fig)
    return str(path)


def sequence_markup(result, max_bases=1600):
    """Standalone HTML mapping primers and mutation on the reference-plus strand."""
    pair = result["selected_primer"]
    if not pair or not result["wt_amplicon"]:
        return "<p>No amplicon available.</p>"
    dna = result["wt_amplicon"]
    full_start = pair["amplicon_start_1based"]
    variant = result["position"] - full_start
    left = len(pair["forward"])
    right_start = len(dna) - len(pair["reverse"])
    selected_start = max(0, variant - max_bases//2) if len(dna) > max_bases else 0
    selected_end = min(len(dna), selected_start + max_bases)
    selected_start = max(0, selected_end - max_bases)
    rows=[]
    for offset in range(selected_start, selected_end, 60):
        fragments=[]
        for i in range(offset,min(offset+60,selected_end)):
            base=escape(dna[i])
            if i == variant:
                kind="variant"
            elif 0 <= i < left:
                kind="forward"
            elif right_start <= i < len(dna):
                kind="reverse"
            else:
                kind="normal"
            fragments.append(f'<span class="{kind}">{base}</span>')
        rows.append(f'<div><span class="loc">{full_start + offset:>9,}</span> ' +
                    "".join(fragments) + "</div>")
    html = """<!DOCTYPE html><html><head><meta charset="utf-8"><style>
body{font-family:system-ui,Arial,sans-serif;color:#172033;margin:18px}
pre{font-family:ui-monospace,Menlo,monospace;font-size:14px;line-height:1.8;
overflow-x:auto;padding:14px;background:#f4f7fc;border:1px solid #e1e8f2;border-radius:10px}
.loc{color:#5e6d83;padding-right:12px}
.forward{background:#b6ebcc;color:#065f46}
.reverse{background:#bfd6ff;color:#1e40af}
.variant{background:#ffe3a0;color:#7c2d12;font-weight:800}
.legend{font-size:13px;margin:12px 0}.tag{padding:4px 8px;margin-right:8px}
</style></head><body>"""
    html += "<h3>Reference-strand amplicon map</h3>"
    html += (f"<p>Input sequence coordinates (1-based): {full_start:,}–"
             f"{pair['amplicon_end_1based']:,}; variant {result['mutation']} "
             f"(position {result['position']:,}).</p>")
    html += ('<div class="legend"><span class="tag forward">Forward primer binding</span>'
             '<span class="tag variant">Variant base (WT)</span>'
             '<span class="tag reverse">Reverse primer binding (plus-strand reverse-complement target)</span></div>')
    html += "<pre>" + "\n".join(rows) + "</pre>"
    html += (f"<p><b>Forward (5′→3′):</b> {escape(pair['forward'])} | "
             f"<b>Reverse (5′→3′):</b> {escape(pair['reverse'])} | "
             f"<b>Product:</b> {pair['product_size']} bp</p>")
    if selected_start or selected_end < len(dna):
        html += "<p>Long amplicon truncated to variant-centered preview; full amplicon is in analysis.json.</p>"
    html += "<p>Primer binding shown schematically; genomic specificity and amplification are not validated.</p></body></html>"
    return html


def render_sequence_figure(result, output_path):
    pair = result["selected_primer"]
    if not pair:
        raise ValueError("Cannot draw primer map without selected primers.")
    begin = pair["amplicon_start_1based"]
    end = pair["amplicon_end_1based"]
    variant = result["position"]
    forward_end = begin + len(pair["forward"]) - 1
    reverse_begin = end - len(pair["reverse"]) + 1
    fig, ax=plt.subplots(figsize=(10, 2.6))
    ax.plot([begin,end],[0,0],color="#334155",lw=6,solid_capstyle="round")
    ax.plot([begin,forward_end],[0,0],color="#16a34a",lw=9)
    ax.plot([reverse_begin,end],[0,0],color="#2563eb",lw=9)
    ax.scatter([variant],[0],color="#f59e0b",s=155,marker="D",zorder=6)
    ax.text(variant,.19,f"{result['mutation']}\n{variant:,} bp",ha="center",fontsize=9)
    ax.text(begin,-.24,f"F 5'→3'\n{pair['forward']}",ha="left",fontsize=8,color="#166534")
    ax.text(end,-.24,f"R 5'→3'\n{pair['reverse']}",ha="right",fontsize=8,color="#1d4ed8")
    ax.set_xlim(begin-8,end+8); ax.set_ylim(-.42,.48); ax.axis("off")
    ax.set_title(f"Primer binding positions | {pair['product_size']} bp amplicon",fontsize=12)
    fig.tight_layout()
    path=Path(output_path);path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=190,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return str(path)
