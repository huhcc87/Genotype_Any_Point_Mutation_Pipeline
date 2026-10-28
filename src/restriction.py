from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, CommOnly, Analysis
from Bio.Restriction import AllEnzymes
import yaml
from pathlib import Path

def load_preferred_enzymes(config_path="config/preferred_enzymes.yaml"):
    """Load preferred enzymes from YAML config."""
    with open(config_path) as f:
        config = yaml.safe_load(f)
    preferred_names = {e["name"] for e in config.get("preferred", [])}
    return preferred_names, config.get("scoring", {})

def find_differential_enzymes(wt_seq: str, mut_seq: str, max_results=15):
    """
    Find restriction enzymes that cut differently between wild-type and mutant.
    Returns a ranked list of useful enzymes.
    """
    preferred_names, scoring = load_preferred_enzymes()
    
    wt = Seq(wt_seq.upper())
    mut = Seq(mut_seq.upper())
    
    # Use commercially available enzymes
    batch = RestrictionBatch(CommOnly)
    
    ana_wt = Analysis(batch, wt, linear=True)
    ana_mut = Analysis(batch, mut, linear=True)
    
    results = []
    
    for enzyme in batch:
        try:
            cuts_wt = sorted(ana_wt.full().get(enzyme, []))
            cuts_mut = sorted(ana_mut.full().get(enzyme, []))
            
            if cuts_wt == cuts_mut:
                continue  # No difference
            
            # Calculate a simple score
            score = 0
            name = str(enzyme)
            
            # Preferred enzyme bonus
            if name in preferred_names:
                score += scoring.get("preferred_bonus", 30)
            
            # Prefer enzymes that cut only once in one allele
            if len(cuts_wt) == 1 and len(cuts_mut) == 0:
                score += scoring.get("single_cut_bonus", 20)
            elif len(cuts_mut) == 1 and len(cuts_wt) == 0:
                score += scoring.get("single_cut_bonus", 20)
            
            # Commercial availability
            score += scoring.get("commercial_bonus", 10)
            
            results.append({
                "enzyme": name,
                "site": str(enzyme.site) if hasattr(enzyme, "site") else "",
                "cuts_wt": cuts_wt,
                "cuts_mut": cuts_mut,
                "score": score,
                "preferred": name in preferred_names
            })
        except Exception:
            continue
    
    # Sort by score (highest first)
    results = sorted(results, key=lambda x: x["score"], reverse=True)
    return results[:max_results]

def print_enzyme_table(results):
    """Pretty print the enzyme results."""
    print("\n{:<12} {:<12} {:<20} {:<20} {:<8}".format(
        "Enzyme", "Site", "Cuts WT", "Cuts Mut", "Score"))
    print("-" * 75)
    for r in results:
        print("{:<12} {:<12} {:<20} {:<20} {:<8}".format(
            r["enzyme"],
            r["site"],
            str(r["cuts_wt"]),
            str(r["cuts_mut"]),
            r["score"]
        ))
