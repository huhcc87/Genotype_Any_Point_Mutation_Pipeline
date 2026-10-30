"""PCR-amplicon-specific, allele-discriminating restriction enzyme screening."""
from pathlib import Path
import yaml
from Bio.Seq import Seq
from Bio.Restriction import CommOnly


def load_preferred_enzymes(config_path=None):
    if config_path is None:
        config_path = Path(__file__).resolve().parents[1] / "config" / "preferred_enzymes.yaml"
    try:
        with open(config_path, encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
    except FileNotFoundError:
        config = {}
    return {x["name"] for x in config.get("preferred", [])}, config.get("scoring", {})


def digest_fragments(amplicon_size: int, cuts: list[int]) -> list[int]:
    """Convert 0-based boundaries between bases to fragment sizes in bp."""
    if amplicon_size <= 0:
        raise ValueError("Amplicon size must be positive.")
    boundaries = [0] + sorted({c for c in cuts if 0 < c < amplicon_size}) + [amplicon_size]
    return sorted((b - a for a, b in zip(boundaries, boundaries[1:])), reverse=True)


def restriction_cuts(enzyme, sequence: str) -> list[int]:
    """Biopython restriction search coordinates are 1-based; convert to 0-based cut boundaries."""
    return sorted({
        site - 1 for site in enzyme.search(Seq(sequence), linear=True)
        if 0 < site - 1 < len(sequence)
    })


def find_differential_enzymes(wt_seq: str, mut_seq: str, max_results=15):
    """Find differentially cleaving enzymes on SAME-LENGTH PCR amplicons, not genomic template.

    This is in-silico cleavage only; methylation, partial digestion and
    commercial enzyme stock/conditions require external validation.
    """
    wt_seq, mut_seq = wt_seq.upper(), mut_seq.upper()
    if not wt_seq or len(wt_seq) != len(mut_seq):
        raise ValueError("WT and mutant amplicons must have equal positive lengths.")
    preferred, scoring = load_preferred_enzymes()
    results = []
    for enzyme in CommOnly:
        try:
            cw = restriction_cuts(enzyme, wt_seq)
            cm = restriction_cuts(enzyme, mut_seq)
            fw = digest_fragments(len(wt_seq), cw)
            fm = digest_fragments(len(mut_seq), cm)
            if fw == fm:
                continue
            # Fragments smaller than 20 bp may not resolve on ordinary agarose gels.
            distinct = sorted(set(fw).symmetric_difference(set(fm)))
            informative = [x for x in distinct if x >= 20]
            score = (scoring.get("preferred_bonus", 30) if str(enzyme) in preferred else 0)
            if (len(cw) == 1 and not cm) or (len(cm) == 1 and not cw):
                score += scoring.get("single_cut_bonus", 20)
            score += scoring.get("commercial_bonus", 10)
            if informative:
                score += scoring.get("good_fragment_separation", 15)
            if not informative:
                score -= 50
            if any(x < 20 for x in fw + fm):
                score -= 10
            results.append({
                "enzyme": str(enzyme),
                "site": str(enzyme.site),
                "cuts_wt": cw,
                "cuts_mut": cm,
                "fragments_wt": fw,
                "fragments_mut": fm,
                "fragments_het": sorted(fw + fm, reverse=True),
                "score": score,
                "preferred": str(enzyme) in preferred,
                "resolvable": bool(informative),
            })
        except (ValueError, TypeError):
            continue
    results.sort(key=lambda x: (x["resolvable"], x["score"], -len(x["cuts_wt"] + x["cuts_mut"])), reverse=True)
    return results[:max_results]


def print_enzyme_table(results):
    print(f"{'Enzyme':<14} {'Site':<12} {'WT fragments':<25} {'Mut fragments':<25} {'Score':<6}")
    print("-" * 90)
    for row in results:
        print(f"{row['enzyme']:<14} {row['site']:<12} {str(row['fragments_wt']):<25} {str(row['fragments_mut']):<25} {row['score']:<6}")
