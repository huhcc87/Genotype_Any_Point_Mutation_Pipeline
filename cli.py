#!/usr/bin/env python3
"""CLI sharing the Streamlit assay engine."""
import csv
import json
from pathlib import Path
import click
from src.assay import analyze_assay
from src.sequence import read_sequence_file
from src.gel import draw_virtual_gel

def save(result,output_dir):
    directory=Path(output_dir)
    directory.mkdir(parents=True,exist_ok=True)
    stem="".join(c if c.isalnum() or c in "-_" else "_" for c in result["name"])[:60]
    stamp=result["timestamp_utc"].replace(":","").replace("-","")
    base=directory / (stem+"_"+stamp)
    target=base.with_suffix(".json")
    target.write_text(json.dumps(result,indent=2),encoding="utf-8")
    if result["enzymes"] and result["selected_primer"]:
        e=result["enzymes"][0]
        draw_virtual_gel(result["selected_primer"]["product_size"],e["cuts_wt"],e["cuts_mut"],e["enzyme"],str(base)+"_gel.png")
    return target

@click.group()
def cli():
    """Research-use PCR-RFLP tools."""

@cli.command()
@click.option("--seq",required=True)
@click.option("--mutation",required=True)
@click.option("--name",default="Variant")
@click.option("--max-enzymes",default=10,type=click.IntRange(1,100))
@click.option("--output-dir",default="output")
def analyze(seq,mutation,name,max_enzymes,output_dir):
    try:
        result=analyze_assay(read_sequence_file(seq) if __import__('pathlib').Path(seq).is_file() else seq,mutation,name,max_enzymes)
        path=save(result,output_dir)
    except (OSError,ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"{result['status']}: {result['mutation']} | {len(result['enzymes'])} enzymes")
    click.echo(str(path))

@cli.command()
@click.option("--csv","csv_path",required=True,type=click.Path(exists=True))
@click.option("--output-dir",default="output")
def batch(csv_path,output_dir):
    errors=0
    with open(csv_path,newline="",encoding="utf-8") as source:
        reader=csv.DictReader(source)
        if not {"name","seq","mutation"}.issubset(reader.fieldnames or []):
            raise click.ClickException("Required CSV columns: name,seq,mutation")
        for i,row in enumerate(reader,1):
            try:
                result=analyze_assay(read_sequence_file(row["seq"]) if __import__('pathlib').Path(row["seq"]).is_file() else row["seq"],row["mutation"],row["name"])
                save(result,output_dir)
                click.echo(f"{i}: {result['status']}")
            except (OSError,ValueError) as exc:
                errors+=1
                click.echo(f"{i}: ERROR {exc}",err=True)
    if errors:
        raise click.ClickException(f"{errors} rows failed")
if __name__=="__main__":
    cli()
