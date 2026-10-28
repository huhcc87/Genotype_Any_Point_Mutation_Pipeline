from pathlib import Path

def generate_html_report(
    name: str,
    mutation_info: dict,
    enzymes: list,
    primers: list,
    gel_path: str = None,
    output_path: str = "report.html"
):
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>PCR-RFLP Report: {name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        img {{ max-width: 600px; }}
    </style>
</head>
<body>
    <h1>PCR-RFLP Design Report: {name}</h1>
    <p><strong>Date:</strong> {mutation_info.get('date', '')}</p>
    
    <h2>Mutation</h2>
    <p>Position: {mutation_info.get('position', 'N/A')} | 
       Change: {mutation_info.get('ref', '?')} → {mutation_info.get('alt', '?')} | 
       Length: {mutation_info.get('seq_length', 'N/A')} bp</p>
"""

    if enzymes:
        best = enzymes[0]
        html += f"<h2>Recommended Enzyme</h2><p>{best['enzyme']} (site: {best['site']})</p>"

    html += "<h2>Primers</h2><ul>"
    for p in primers[:3]:
        html += f"<li>Pair {p['pair']}: {p['forward']} / {p['reverse']} ({p['product_size']} bp)</li>"
    html += "</ul>"

    if gel_path and Path(gel_path).exists():
        html += f'<h2>Virtual Gel</h2><img src="{gel_path}" alt="Virtual Gel">'

    html += "</body></html>"

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)
    return output_path
