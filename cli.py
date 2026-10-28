#!/usr/bin/env python3
"""
Genotype Pipeline CLI
Design PCR-RFLP assays for point mutations.
"""

import click
from pathlib import Path
from datetime import datetime

from src.sequence import load_sequence, parse_mutation, apply_mutation
from src.restriction import find_differential_enzymes, print_enzyme_table
from src.primers import design_primers, simple_end_primers
from src.gel import draw_virtual_gel
from src.report import generate_markdown_report


@click.command()
@click.option("--seq", required=True, help="DNA sequence or path to FASTA file")
@click.option("--mutation", required=True, help="Mutation e.g. G1237A or c.G1237A")
@click.option("--name", default="Variant", help="Name of the variant")
@click.option("--max-enzymes", default=5, type=int, help="Maximum number of enzymes to show")
def analyze(seq, mutation, name, max_enzymes):
    """Design PCR-RFLP assay for a point mutation."""
    click.echo("=" * 60)
    click.echo("  Genotype Pipeline - PCR-RFLP Design")
    click.echo("=" * 60)

    # 1. Load sequence
    try:
        sequence = load_sequence(seq)
        click.echo(f"\n[1] Sequence loaded: {len(sequence)} bp")
    except Exception as e:
        click.echo(f"Error loading sequence: {e}")
        return

    # 2. Parse mutation
    try:
        pos, ref, alt = parse_mutation(mutation)
        click.echo(f"[2] Mutation parsed: position {pos+1} ({ref}>{alt})")
    except Exception as e:
        click.echo(f"Error parsing mutation: {e}")
        return

    # 3. Apply mutation
    try:
        wt_seq, mut_seq = apply_mutation(sequence, pos, ref, alt)
        click.echo("[3] Wild-type and mutant sequences created")
        ctx_start = max(0, pos - 8)
        ctx_end = min(len(sequence), pos + 9)
        click.echo(f"    Context: ...{wt_seq[ctx_start:ctx_end]}...")
        click.echo(f"    Mutant : ...{mut_seq[ctx_start:ctx_end]}...")
    except Exception as e:
        click.echo(f"Error applying mutation: {e}")
        return

    # 4. Find enzymes
    click.echo("\n[4] Searching for differential restriction enzymes...")
    enzymes = find_differential_enzymes(wt_seq, mut_seq, max_results=max_enzymes)
    if enzymes:
        print_enzyme_table(enzymes)
        best = enzymes[0]
        click.echo(f"\n  Best enzyme: {best['enzyme']} (site: {best['site']})")
    else:
        click.echo("  No differential enzymes found.")

    # 5. Design primers
    click.echo("\n[5] Designing primers...")
    primers = design_primers(sequence, pos)
    if not primers:
        click.echo("  Primer3 found no pairs — using simple end primers...")
        primers = simple_end_primers(sequence)

    if primers:
        for p in primers[:3]:
            click.echo(f"\n  Pair {p['pair']}:")
            click.echo(f"    Forward : {p['forward']} (Tm={p['forward_tm']}°C)")
            click.echo(f"    Reverse : {p['reverse']} (Tm={p['reverse_tm']}°C)")
            click.echo(f"    Product : {p['product_size']} bp")

    # 6. Virtual gel
    gel_path = None
    if enzymes and primers:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_name = name.replace(" ", "_")
        gel_file = Path("output") / f"{safe_name}_{timestamp}_gel.png"
        try:
            gel_path = draw_virtual_gel(
                product_size=primers[0]["product_size"],
                enzyme=enzymes[0],["enzyme"],
                output_path=str(gel_file)
            )
            click.echo(f"\n[6] Virtual gel saved: {gel_path}")
        except Exception as e:
            click.echo(f"  Could not generate gel: {e}")

    # 7. Report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = name.replace(" ", "_")
    report_file = Path("output") / f"{safe_name}_{timestamp}_report.md"
    try:
        generate_markdown_report(
            name=name,
            mutation_info={
                "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
                "position": pos + 1,
                "ref": ref,
                "alt": alt,
                "seq_length": len(sequence),
            },
            enzymes=enzymes,
            primers=primers,
            gel_path=str(gel_path) if gel_path else None,
            output_path=str(report_file),
        )
        click.echo(f"[7] Report saved: {report_file}")
    except Exception as e:
        click.echo(f"Could not generate report: {e}")

    click.echo("\n" + "=" * 60)
    click.echo("  Analysis complete")
    click.echo("=" * 60)


if __name__ == "__main__":
    analyze()
