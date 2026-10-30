"""Heuristic ranking of computed PCR-RFLP patterns. Not empirical enzyme performance."""
def rank_enzymes(enzymes):
    ranked = []
    for e in enzymes:
        wt = e.get("fragments_wt", [])
        mut = e.get("fragments_mut", [])
        if not wt or not mut:
            continue
        unique_wt = [x for x in wt if x not in mut]
        unique_mut = [x for x in mut if x not in wt]
        diff_visible = [x for x in unique_wt + unique_mut if x >= 50]
        score = 0
        reasons = []
        if unique_wt and unique_mut:
            score += 30
            reasons.append("alleles have distinct predicted fragments")
        elif unique_wt or unique_mut:
            score += 20
            reasons.append("one allele has a distinguishing fragment")
        if diff_visible:
            score += 25
            reasons.append("distinguishing fragment at least 50 bp")
        else:
            reasons.append("differential fragments below practical visualization threshold")
        if len(e.get("cuts_wt", [])) <= 1 and len(e.get("cuts_mut", [])) <= 1:
            score += 20
            reasons.append("simple restriction pattern")
        else:
            reasons.append("multiple restriction sites")
        if all(x >= 50 for x in wt + mut):
            score += 15
        else:
            reasons.append("at least one fragment below 50 bp")
        if e.get("preferred"):
            score += 10
            reasons.append("user-configured preferred enzyme")
        # Distinct bands must be separated from all other genotype bands.
        band_sizes = sorted(set(wt + mut))
        if diff_visible and any(
            any(0 < abs(x - y) < max(15, x * .08) for y in band_sizes if y != x)
            for x in diff_visible
        ):
            score -= 15
            reasons.append("closely spaced fragments may co-migrate")
        score = max(0, min(100, score))
        if score >= 80:
            grade = "A"
        elif score >= 65:
            grade = "B"
        elif score >= 45:
            grade = "C"
        else:
            grade = "D"
        entry = dict(e)
        entry.update({
            "priority_score": score,
            "design_grade": grade,
            "design_notes": "; ".join(reasons),
            "enzyme_status": "PREDICTED ONLY",
        })
        ranked.append(entry)
    ranked.sort(key=lambda e: (e["priority_score"], -len(e["fragments_wt"] + e["fragments_mut"]), e["enzyme"]), reverse=True)
    for i, e in enumerate(ranked, 1):
        e["priority_rank"] = i
    return ranked
