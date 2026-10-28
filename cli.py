#!/usr/bin/env python3
"""
Genotype Pipeline - Reusable PCR-RFLP design tool
https://github.com/yourusername/Genotype_Pipeline
"""

import click
from pathlib import Path
import sys
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent / "src"))

from sequence import load_sequence, apply_mutation, parse_mutation
from restriction import find_differential_enzymes, print_enzyme_table
from primers import design_primers, simple_end_primers
from gel import draw_virtual_gel
from report import generate_markdown_report

@click.command()
@click.option("--seq", required=True, help="DNA sequence or path to FASTA file")
@click.option("--mutation", required=True, help="Mutation e.g. G1237A or c.G1237A")
@click.option("--name", default="Variant", help="Name of the variant")
@click.option("--max-enzymes", default=10, help="Maximum enzymes to show")
@click.option("--output-dir", default="output", help="Output directory")
@click.option("--no-gel", is_flag=True, help="Skip virtual gel generation")
def analyze(seq, mutation, name, max_enzymes, output_dir, no_gel):
    """Design PCR-RFLP assay for a point mutation."""

    click.echo("\n" + "="*60)
    click.echo("  Genotype Pipeline - PCR-RFLP Design")
    click.echo("="*60)

    # 1. Load sequence
    try:
        sequence = load_sequence(seq)
        click.echo(f"\n[1] Sequence loaded: {len(sequence)} bp")
    except Exception as e:
        click.echo(f"Error loading sequence: {e}", err=True)
        return

    # 2. Parse mutation
    try:
        pos, ref, alt = parse_mutation(mutation)
        click.echo(f"[2] Mutation parsed: position {pos+1} ({ref}>{alt})")
    except Exception as e:
        click.echo(f"Error parsing mutation: {e}", err=True)
        return

    # 3. Apply mutation
    try:
        wt, mut = apply_mutation(sequence, pos, ref, alt)
        click.echo(f"[3] Wild-type and mutant sequences created")
        click.echo(f"    Context: ...{wt[max(0,pos-8):pos+10]}...")
        click.echo(f"    Mutant : ...{mut[max(0,pos-8):pos+10]}...")
    except Exception as e:
        click.echo(f"Error applying mutation: {e}", err=True)
        return

    # 4. Find enzymes
    click.echo(f"\n[4] Searching for differential restriction enzymes...")
    enzymes = find_differential_enzymes(wt, mut, max_results=max_enzymes)

    if not enzymes:
        click.echo("  No differential enzymes found.")
        best = None
    else:
        print_enzyme_table(enzymes)
        best = enzymes[0]
        click.echo(f"\n  Best enzyme: {best['enzyme']} (site: {best['site']})")

    # 5. Design primers
    click.echo(f"\n[5] Designing primers...")
    primers = design_primers(sequence, pos)
    if not primers:
        click.echo("  Primer3 found no pairs — using simple end primers...")
        primers = simple_end_primers(sequence)

    if primers:
        for p in primers[:3]:
            click.echo(f"\n  Pair {p['pair']}:")
            click.echo(f"    Forward : {p['forward']}  (Tm={p['forward_tm']}°C)")
            click.echo(f"    Reverse : {p['reverse']}  (Tm={p['reverse_tm']}°C)")
            click.echo(f"    Product : {p['product_size']} bp")
    else:
        click.echo("  Could not design primers.")

    # 6. Virtual gel + Report
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = name.replace(" ", "_")

    gel_path = None
    if best and not no_gel and primers:
        amplicon = primers[0]["product_size"]
        gel_file = output_path / f"{safe_name}_{timestamp}_gel.png"
        try:
            gel_path = draw_virtual_gel(
                amplicon_size=amplicon,
                cuts_wt=best["cuts_wt"],
                cuts_mut=best["cuts_mut"],
                enzyme_name=best["enzyme"],
                output_path=str(gel_file),
                title=f"{name} – {best['enzyme']}"
            )
            click.echo(f"\n[6] Virtual gel saved: {gel_path}")
        except Exception as e:
            click.echo(f"  Could not generate gel: {e}")

    # Report
    report_file = output_path / f"{safe_name}_{timestamp}_report.md"
    try:
        generate_markdown_report(
            name=name,
            sequence=sequence,
            mutation=mutation,
            pos=pos,
            ref=ref,
            alt=alt,
            enzymes=enzymes,
            primers=primers,
            gel_path=str(gel_path) if gel_path else None,
            output_path=str(report_file)
        )
        click.echo(f"[7] Report saved: {report_file}")
    except Exception as e:
        click.echo(f"  Could not generate report: {e}")

    click.echo("\n" + "="*60)
    click.echo("  Analysis complete")
    click.echo("="*60 + "\n")

if __name__ == "__main__":
    analyze()
