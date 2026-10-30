"""DNA ladder presets. Vendor-specific sizes cited to manufacturer sources.

All simulated migration and brightness remain illustrative, not measured.
"""
LADDERS = {
    "Thermo GeneRuler 100 bp (SM0241)": {
        "vendor": "Thermo Fisher Scientific", "catalog":"SM0241",
        "sizes": [1000,900,800,700,600,500,400,300,200,100],
        "source": "https://tools.thermofisher.com/content/sfs/manuals/MAN0012995_GeneRuler_100bp_DNALadder_50ug_UG.pdf",
        "note":"Nominal marker sizes; 500 bp reference band. Not calibrated intensity."
    },
    "Thermo GeneRuler 50 bp (SM0373)": {
        "vendor": "Thermo Fisher Scientific", "catalog":"SM0373",
        "sizes": [1000,900,800,700,600,500,400,300,250,200,150,100,50],
        "source": "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013018_GeneRuler_50bp_DNALadder_RTU_50ug_UG.pdf",
        "note":"Manufacturer-reported nominal sizes; 500 and 250 bp orientation bands."
    },
    "NEB Quick-Load Purple 100 bp (N0551)": {
        "vendor":"New England Biolabs","catalog":"N0551",
        "sizes":[1517,1200,1000,900,800,700,600,517,500,400,300,200,100],
        "source":"https://www.neb.com/en-us/products/n0551-quick-load-purple-100-bp-dna-ladder",
        "note":"500 and 517 bp marker bands may be indistinguishable on ordinary agarose."
    },
    "NEB 1 kb DNA Ladder (N3232)": {
        "vendor":"New England Biolabs", "catalog":"N3232",
        "sizes":[10002,8001,6001,5001,4001,3001,2000,1500,1000,517,500],
        "source":"https://www.neb.com/en-us/products/n3232-1-kb-dna-ladder",
        "note":"Manufacturer-reported fragment sizes. Designed for 0.5–10 kb, not for ~50 bp PCR-RFLP fragments."
    },
    "Generic 25 bp (illustration)": {
        "vendor":"Generic, not vendor-specific","catalog":"Illustration",
        "sizes":[500,450,400,350,300,275,250,225,200,175,150,125,100,75,50,25],
        "source":"",
        "note":"Illustrative reference intervals only; not a real product specification."
    },
    "Generic 50 bp (illustration)": {
        "vendor":"Generic, not vendor-specific","catalog":"Illustration",
        "sizes":[1000,900,800,700,600,500,450,400,350,300,250,200,150,100,50],
        "source":"",
        "note":"Illustrative reference intervals only; not a real product specification."
    },
    "Generic 100 bp (illustration)": {
        "vendor":"Generic, not vendor-specific","catalog":"Illustration",
        "sizes":[1500,1200,1000,900,800,700,600,500,400,300,200,100],
        "source":"",
        "note":"Illustrative reference intervals only; not a real product specification."
    },
    "Generic 1 kb (illustration)": {
        "vendor":"Generic, not vendor-specific","catalog":"Illustration",
        "sizes":[10000,8000,6000,5000,4000,3000,2000,1500,1000,750,500,250],
        "source":"",
        "note":"Illustrative reference intervals only; not a real product specification."
    },
}

def ladder_notes(ladder_key, fragments):
    ladder=LADDERS[ladder_key]
    lo,hi=min(ladder["sizes"]),max(ladder["sizes"])
    outside=sorted({int(x) for x in fragments if x < lo or x > hi})
    note=ladder["note"]
    if outside:
        note+=" Predicted sample fragments outside marker range: "+", ".join(map(str,outside))+" bp."
    if any(x<50 for x in fragments):
        note+=" Fragments smaller than 50 bp may be difficult to resolve on standard agarose."
    return note

def recommended_ladders(fragments):
    """Rank by coverage and distance of fragments to marker range; no performance claim."""
    if not fragments:
        return list(LADDERS)
    low=min(fragments);high=max(fragments)
    def sortkey(item):
        label,record=item
        lo,hi=min(record["sizes"]),max(record["sizes"])
        coverage=sum(lo<=x<=hi for x in fragments)
        overhang=max(0,lo-low)+max(0,high-hi)
        return (-coverage,overhang,hi-lo,label)
    return [name for name,_ in sorted(LADDERS.items(),key=sortkey)]
