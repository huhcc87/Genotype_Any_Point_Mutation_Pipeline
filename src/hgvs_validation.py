"""Transcript-aware SNV HGVS validation against the SELECTED coding transcript.

Only coding SNVs c.NRef>Alt are accepted. Protein notation is a check,
not a genomic coordinate; all mappings are validated against reference CDS.
"""
import re
from Bio.Seq import Seq
from .gene_lookup import coding_map

SNV=re.compile(r"^(?:[A-Za-z0-9_.]+:)?c\.(\d+)([ACGT])>([ACGT])$",re.I)
LEGACY=re.compile(r"^c\.([ACGT])(\d+)([ACGT])$",re.I)
PROTEIN=re.compile(r"^(?:p\.)?([A-Z])([1-9]\d*)([A-Z*])$",re.I)

def parse_coding_snv(text):
    value=str(text).strip().replace(" ","")
    match=SNV.fullmatch(value)
    if match: n,ref,alt=int(match[1]),match[2].upper(),match[3].upper()
    else:
        old=LEGACY.fullmatch(value)
        if not old: raise ValueError("Use c.1237G>A (or legacy c.G1237A) for a coding single-base variant.")
        n,ref,alt=int(old[2]),old[1].upper(),old[3].upper()
    if n<1 or ref==alt: raise ValueError("Variant position must be positive and REF must differ from ALT.")
    return n,ref,alt

def validate_coding_snv(transcript, cds, annotation, expected_protein=None):
    pos,ref,alt=parse_coding_snv(annotation)
    if pos>len(cds): raise ValueError(f"c.{pos} exceeds {len(cds)}-bp CDS for selected transcript.")
    if cds[pos-1]!=ref:
        raise ValueError(f"Reference mismatch: submitted c.{pos}{ref}>{alt}, but selected transcript has {cds[pos-1]} at c.{pos}. Select correct accession/version.")
    codon_number=(pos-1)//3+1
    first=(codon_number-1)*3
    wt=cds[first:first+3]
    offset=(pos-1)%3
    mutant=wt[:offset]+alt+wt[offset+1:]
    protein_before=str(Seq(wt).translate())
    protein_after=str(Seq(mutant).translate())
    predicted=f"{protein_before}{codon_number}{protein_after}"
    if expected_protein:
        p=PROTEIN.fullmatch(str(expected_protein).strip().replace(" ",""))
        if not p: raise ValueError("Protein check must use p.E413K or E413K (one-letter amino acids).")
        if (p[1].upper(),int(p[2]),p[3].upper())!=(protein_before,codon_number,protein_after):
            raise ValueError(f"Protein mismatch: expected {expected_protein}, selected transcript predicts p.{predicted}.")
    coords,exons=coding_map(transcript,cds)
    strand=transcript["strand"]
    genomic_ref=str(Seq(ref).complement()) if strand==-1 else ref
    genomic_alt=str(Seq(alt).complement()) if strand==-1 else alt
    return {"status":"REFERENCE_AND_PROTEIN_VERIFIED" if expected_protein else "REFERENCE_VERIFIED",
            "transcript":transcript.get("id"),"coding_hgvs":f"c.{pos}{ref}>{alt}",
            "protein_hgvs":f"p.{predicted}","codon_number":codon_number,
            "reference_codon":wt,"alternate_codon":mutant,
            "chromosome":transcript["seq_region_name"],"genomic_position":coords[pos-1],
            "genomic_ref":genomic_ref,"genomic_alt":genomic_alt,
            "exon_numbers":sorted(set(exons[first:first+3]))}

def validate_batch(transcript,cds,text):
    lines=[ln.strip() for ln in str(text).splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if len(lines)>500: raise ValueError("Maximum 500 coding variants per batch.")
    records=[]
    for line in lines:
        fields=[x.strip() for x in line.replace("\t",",").split(",")]
        try:
            result=validate_coding_snv(transcript,cds,fields[0],fields[1] if len(fields)>1 and fields[1] else None)
            result["input"]=line
        except (ValueError,KeyError,TypeError) as exc:
            result={"input":line,"status":"MISMATCH_OR_UNSUPPORTED","error":str(exc)}
        records.append(result)
    return records
