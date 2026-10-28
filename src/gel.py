import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

def draw_virtual_gel(amplicon_size, cuts_wt, cuts_mut, enzyme_name, 
                     output_path="output/virtual_gel.png", title=None):
    """
    Draw a simple virtual agarose gel showing expected bands.
    """
    fig, ax = plt.subplots(figsize=(6, 8))
    
    # Lanes
    lanes = ["Ladder", "WT undig", "WT dig", "Het dig", "Mut dig"]
    lane_x = [1, 2.5, 4, 5.5, 7]
    
    # Background
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, max(amplicon_size * 1.2, 300))
    ax.set_facecolor("#1a1a2e")
    
    # Draw gel wells
    for x in lane_x:
        ax.add_patch(patches.Rectangle((x-0.3, max(amplicon_size*1.15, 280)), 0.6, 8, 
                                        facecolor="white", edgecolor="gray"))
    
    def draw_band(x, size, intensity=1.0):
        if size <= 0:
            return
        y = size
        height = max(3, min(12, size / 20))
        color = (0.2, 0.9, 0.3, intensity)
        ax.add_patch(patches.Rectangle((x-0.25, y - height/2), 0.5, height,
                                        facecolor=color, edgecolor="none"))
    
    # Ladder (approximate)
    for s in [50, 100, 150, 200, 250, 300, 400, 500]:
        if s < max(amplicon_size * 1.2, 300):
            draw_band(lane_x[0], s, 0.6)
    
    # WT undigested
    draw_band(lane_x[1], amplicon_size)
    
    # WT digested
    if cuts_wt:
        prev = 0
        for c in sorted(cuts_wt) + [amplicon_size]:
            frag = c - prev
            if frag > 5:
                draw_band(lane_x[2], frag)
            prev = c
    else:
        draw_band(lane_x[2], amplicon_size)
    
    # Heterozygote (both patterns)
    draw_band(lane_x[3], amplicon_size, 0.7)
    if cuts_wt:
        prev = 0
        for c in sorted(cuts_wt) + [amplicon_size]:
            frag = c - prev
            if frag > 5:
                draw_band(lane_x[3], frag, 0.7)
            prev = c
    
    # Mutant digested
    if cuts_mut:
        prev = 0
        for c in sorted(cuts_mut) + [amplicon_size]:
            frag = c - prev
            if frag > 5:
                draw_band(lane_x[4], frag)
            prev = c
    else:
        draw_band(lane_x[4], amplicon_size)
    
    # Labels
    ax.set_xticks(lane_x)
    ax.set_xticklabels(lanes, fontsize=9, color="white")
    ax.set_ylabel("Fragment size (bp)", color="white")
    ax.tick_params(colors="white")
    ax.set_title(title or f"Virtual Gel – {enzyme_name}", color="white", fontsize=12)
    
    # Invert y so large fragments are at top
    ax.invert_yaxis()
    
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, facecolor="#1a1a2e")
    plt.close()
    return output_path
