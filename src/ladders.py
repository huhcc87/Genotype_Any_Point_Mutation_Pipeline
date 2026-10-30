"""Verified product-specific or explicitly generic DNA ladder presets.

Vendor metadata links are publicly accessible manufacturer product pages.
Do not treat these as calibrated intensity or migration models.
"""
LADDERS = {
    "Thermo GeneRuler 100 bp (SM0241)": {
        "vendor": "Thermo Fisher Scientific",
        "catalog": "SM0241",
        "sizes": [1000,900,800,700,600,500,400,300,200,100],
        "source": "https://tools.thermofisher.com/content/sfs/manuals/MAN0012995_GeneRuler_100bp_DNALadder_50ug_UG.pdf",
        "note": "500 bp orientation band. Nominal fragment sizes from manufacturer."
    },
    "NEB Quick-Load Purple 100 bp (N0551)": {
        "vendor": "New England Biolabs",
        "catalog": "N0551",
        "sizes": [1517,1200,1000,900,800,700,600,517,500,400,300,200,100],
        "source": "https://www.neb.com/en-us/products/n0551-quick-load-purple-100-bp-dna-ladder",
        "note": "500/517 bp products co-migrate; drawn as separate close sizes only for illustration."
    },
    "Generic 50 bp illustration": {
        "vendor": "Not product-specific",
        "catalog": "Illustrative",
        "sizes": list(range(500,49,-50)),
        "source": "",
        "note": "Illustrative 50 bp intervals, not a manufacturer's validated ladder."
    },
    "Generic 1 kb illustration": {
        "vendor": "Not product-specific",
        "catalog": "Illustrative",
        "sizes": [10000,8000,6000,5000,4000,3000,2000,1500,1000,750,500,250],
        "source": "",
        "note": "Illustrative ladder, not a particular vendor's product."
    },
}


def ladder_notes(ladder_key, fragments):
    ladder = LADDERS[ladder_key]
    lo, hi = min(ladder["sizes"]), max(ladder["sizes"])
    outside = sorted(set(int(x) for x in fragments if x < lo or x > hi))
    note = ladder["note"]
    if outside:
        note += " Some predicted products fall outside nominal marker range: " + ", ".join(map(str, outside)) + " bp."
    if any(x < 50 for x in fragments):
        note += " Predicted fragments <50 bp may be difficult to visualize on standard agarose."
    return note
