"""Illustrative PCR-RFLP agarose gel; bands reflect true diploid allele digests."""
from collections import Counter
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from .restriction import digest_fragments


def gel_fragment_profiles(amplicon_size: int, cuts_wt, cuts_mut):
    wt = digest_fragments(amplicon_size, cuts_wt)
    mut = digest_fragments(amplicon_size, cuts_mut)
    return {
        "WT undigested": [amplicon_size],
        "WT digested": wt,
        "Heterozygote": wt + mut,
        "Mutant digested": mut,
    }


def draw_virtual_gel(amplicon_size, cuts_wt, cuts_mut, enzyme_name,
                     output_path="output/virtual_gel.png", title=None):
    profiles = gel_fragment_profiles(amplicon_size, cuts_wt, cuts_mut)
    ladder = [50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
    ladder = [x for x in ladder if x <= max(amplicon_size * 1.3, 300)]
    labels = ["Ladder"] + list(profiles)
    ymax = max(1000, amplicon_size * 1.5)
    bottom = 20
    def ypos(size):
        return 0.12 + 0.78 * (math.log(ymax) - math.log(max(size, bottom))) / (math.log(ymax) - math.log(bottom))
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("#142033")
    ax.set_facecolor("#142033")
    ax.set_xlim(0.5, len(labels) + 0.5)
    ax.set_ylim(1.03, 0)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=12, ha="right", color="white", fontsize=9)
    ax.set_yticks([ypos(x) for x in (50, 100, 200, 500, 1000) if x <= ymax])
    ax.set_yticklabels([str(x) for x in (50, 100, 200, 500, 1000) if x <= ymax], color="white")
    ax.set_ylabel("Approximate fragment size (bp)", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("#6986a6")
    for lane, sizes in enumerate([ladder] + list(profiles.values()), start=1):
        ax.add_patch(patches.Rectangle((lane - 0.28, 0.065), 0.56, 0.02, color="#8292a5"))
        counts = Counter(sizes)
        for size, copies in counts.items():
            # Short (<20 bp) products are not represented.
            if size < bottom:
                continue
            brightness = min(1.0, 0.40 + 0.16 * copies + 0.35 * (size / max(amplicon_size, 1)))
            ax.add_patch(patches.Rectangle(
                (lane - 0.27, ypos(size) - 0.006), 0.54, 0.012,
                color=(0.25, 0.95, 0.65, brightness)))
    ax.set_title(title or f"Predicted complete digest: {enzyme_name}", color="white")
    fig.tight_layout()
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)
    return str(path)
