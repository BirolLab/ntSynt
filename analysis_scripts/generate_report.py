#!/usr/bin/env python3
"""
generate_report.py

Generate a self-contained HTML report summarising a completed ntSynt
multi-genome synteny run.  All images are embedded as base64 and all
tables are inlined, so the output is a single portable .html file.

Usage (standalone):
    python generate_report.py \
        --abyss-fac        <date>_assemblies/<group>_abyss_fac-summary.tsv \
        --block-stats      ntsynt_run/ntSynt.k24.w1000.synteny_blocks.stats.tsv \
        --ribbon-plot      ntsynt_run/ntsynt-viz/<group>_ribbon-plot.png \
        --discontinuity    ntsynt_run/ntSynt.k24.w1000.discontinuity_reasons.tsv \
        --mash-plot        mash/mash_divergence_boxplot.png \
        --group            lucinidae \
        --output           lucinidae_ntsynt_report.html

Called by the Snakemake rule generate_report via the same arguments.
"""

import argparse
import base64
import csv
import html
import os
import sys
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def encode_image(path: str) -> str:
    """Return a base64-encoded data URI for a PNG image."""
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode()
    return f"data:image/png;base64,{data}"

def format_genome_size(bp):
    """Format the genome size, and return the most appropriate units"""
    if bp >= 1_000_000_000:
        return f"{bp / 1_000_000_000:.1f}", "Gbp", 1_000_000_000
    elif bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f}", "Mbp", 1_000_000
    elif bp >= 1_000:
        return f"{bp / 1_000:.1f}", "kbp", 1_000
    else:
        return str(bp), "bp", 1


def tsv_to_html_table(path: str, table_id: str = "", rename: dict = None) -> str:
    """Read a TSV and return an HTML <table> string.
    
    Args:
        path:     Path to the TSV file.
        table_id: Optional HTML id attribute for the <table> element.
        rename:   Optional dict mapping original column names to display names.
    """
    rename = rename or {}
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    if not rows:
        return "<p><em>No data.</em></p>"

    id_attr = f' id="{table_id}"' if table_id else ""
    lines = [f'<table{id_attr} class="data-table">']

    # Header — apply renames
    lines.append("<thead><tr>")
    for cell in rows[0]:
        display = rename.get(cell, cell)
        lines.append(f"  <th>{html.escape(display)}</th>")
    lines.append("</tr></thead>")

    # Body — unchanged
    lines.append("<tbody>")
    for row in rows[1:]:
        lines.append("<tr>")
        for cell in row:
            lines.append(f"  <td>{html.escape(cell)}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")

    return "\n".join(lines)


def section(title: str, content: str, section_id: str = "") -> str:
    id_attr = f' id="{section_id}"' if section_id else ""
    return f"""
    <section{id_attr} class="report-section">
      <h2>{title}</h2>
      {content}
    </section>"""


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ntSynt Report — {group}</title>
  <style>
    /* ------------------------------------------------------------------ */
    /* Design: scientific report — dark slate + teal accent, monospace     */
    /* tables, clean typographic hierarchy                                  */
    /* ------------------------------------------------------------------ */

    @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500&family=Inter:wght@300;400;600&display=swap');

    :root {{
      --bg:          #f8fafc;
      --surface:     #ffffff;
      --surface2:    #f1f5f9;
      --border:      #e2e8f0;
      --accent:      #0f6eb4;
      --accent-dim:  #dbeafe;
      --text:        #1e293b;
      --text-muted:  #64748b;
      --heading:     #0f172a;
      --danger:      #dc2626;
      --mono:        'JetBrains Mono', monospace;
      --sans:        'Inter', sans-serif;
      --radius:      8px;
    }}

    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

    html {{ scroll-behavior: smooth; }}

    body {{
      background: var(--bg);
      color: var(--text);
      font-family: var(--sans);
      font-size: 15px;
      line-height: 1.65;
    }}

    /* --- header -------------------------------------------------------- */
    header {{
      background: var(--surface);
      border-bottom: 1px solid var(--border);
      padding: 2.5rem 3rem;
      display: flex;
      align-items: flex-end;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 1rem;
    }}

    header h1 {{
      font-family: var(--mono);
      font-size: 2rem;
      font-weight: 500;
      color: var(--heading);
      letter-spacing: -0.02em;
    }}

    header h1 span {{
      color: var(--accent);
    }}

    .meta {{
      font-family: var(--mono);
      font-size: 0.78rem;
      color: var(--text-muted);
      text-align: right;
      line-height: 1.8;
    }}

    /* --- nav ----------------------------------------------------------- */
    nav {{
      background: var(--surface2);
      border-bottom: 1px solid var(--border);
      padding: 0 3rem;
      display: flex;
      gap: 0;
      overflow-x: auto;
    }}

    nav a {{
      display: block;
      padding: 0.85rem 1.2rem;
      font-family: var(--mono);
      font-size: 1rem;
      color: var(--text-muted);
      text-decoration: none;
      border-bottom: 2px solid transparent;
      white-space: nowrap;
      transition: color 0.15s, border-color 0.15s;
    }}

    nav a:hover {{
      color: var(--accent);
      border-bottom-color: var(--accent);
    }}

    /* --- layout -------------------------------------------------------- */
    main {{
      max-width: 1100px;
      margin: 0 auto;
      padding: 2.5rem 3rem 5rem;
    }}

    .report-section {{
      margin-bottom: 3.5rem;
      padding-top: 1rem;
    }}

    .report-section h2 {{
      font-family: var(--mono);
      font-size: 1.2rem;
      font-weight: 500;
      letter-spacing: 0.12em;
      text-transform: uppercase;
      color: var(--accent);
      border-left: 3px solid var(--accent);
      padding-left: 0.75rem;
      margin-bottom: 1.5rem;
    }}

    /* --- tables -------------------------------------------------------- */
    .table-wrap {{
      overflow-x: auto;
      border: 1px solid var(--border);
      border-radius: var(--radius);
    }}

    table.data-table {{
      width: 100%;
      border-collapse: collapse;
      font-family: var(--mono);
      font-size: 0.9rem;
    }}

    table.data-table thead {{
      background: var(--surface2);
    }}

    table.data-table th {{
      padding: 0.65rem 1rem;
      text-align: left;
      color: var(--accent);
      font-weight: 500;
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
    }}

    table.data-table td {{
      padding: 0.55rem 1rem;
      border-bottom: 1px solid var(--border);
      color: var(--text);
    }}

    table.data-table tr:last-child td {{
      border-bottom: none;
    }}

    table.data-table tbody tr:hover {{
      background: var(--surface2);
    }}

    .table-note {{
      margin: 0.75rem 0 0 0;
      font-size: 0.7rem;
      color: var(--text-muted);
      font-family: var(--mono);
    }}

    /* --- figures ------------------------------------------------------- */
    .figure-wrap {{
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      padding: 1.5rem;
      display: inline-block;
      max-width: 100%;
    }}

    .figure-wrap img {{
      display: block;
      max-width: 100%;
      height: auto;
      border-radius: 4px;
    }}

    figcaption {{
      margin-top: 0.75rem;
      font-size: 0.78rem;
      color: var(--text-muted);
      font-family: var(--mono);
    }}

    /* --- two-column figures ------------------------------------------- */
    .figures-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
      gap: 1.5rem;
    }}

    /* --- footer -------------------------------------------------------- */
    footer {{
      background: var(--surface);
      border-top: 1px solid var(--border);
      padding: 1.5rem 3rem;
      font-family: var(--mono);
      font-size: 0.75rem;
      color: var(--text-muted);
      text-align: center;
    }}

  </style>
</head>
<body>

<header>
  <h1>ntSynt report &mdash; <span>{group}</span></h1>
  <div class="meta">
    <div>Generated {date}</div>
    <div>ntSynt multi-genome synteny pipeline</div>
  </div>
</header>

<nav>
  <a href="#assembly-stats">Assembly stats</a>
  <a href="#synteny-blocks">Synteny blocks</a>
  <a href="#discontinuity">Discontinuity reasons</a>
  <a href="#divergence">Mash divergence</a>
  <a href="#ribbon">Ribbon plot</a>
</nav>

<main>
  {sections}
</main>

<footer>
  ntSynt pipeline report &bull; {group} &bull; {date}
</footer>

</body>
</html>
"""


# ---------------------------------------------------------------------------
# Build report
# ---------------------------------------------------------------------------

def build_report(args: argparse.Namespace) -> str:
    group_display = args.group.replace("_", " ").capitalize()
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M")

    sections_html = []

# 1. Assembly stats (abyss-fac)
    if args.abyss_fac and os.path.exists(args.abyss_fac):
        with open(args.abyss_fac, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader)  # single summary row

        n_count   = row["n_count"]
        chr_range = f"{row['n_min']} – {row['n_max']}"
        size_min_val, size_unit, min_divisor = format_genome_size(int(row["sum_min"]))

        size_range = f"{size_min_val} – {int(row['sum_max']) / min_divisor:.1f}"

        summary_rows = [
            ("Number of assemblies", n_count),
            ("Range of number of chromosomes", chr_range),
            (f"Genome size range ({size_unit})", size_range),
        ]

        table_lines = [
            '<table id="tbl-abyss" class="data-table">',
            "<thead><tr><th>Statistic</th><th>Value</th></tr></thead>",
            "<tbody>",
        ]
        for label, value in summary_rows:
            table_lines.append(
                f"<tr><td>{html.escape(label)}</td><td>{html.escape(str(value))}</td></tr>"
            )
        table_lines.append("</tbody></table>")

        content = (
            f'<div class="table-wrap">{"".join(table_lines)}</div>'
            f'<p class="table-note">Source: {html.escape(args.abyss_fac)}</p>'
        )
    else:
        content = "<p><em>abyss-fac summary not found.</em></p>"
    sections_html.append(section("Assembly statistics", content, "assembly-stats"))
    
# 2. Synteny block stats
    if args.block_stats and os.path.exists(args.block_stats):
        BLOCK_STATS_COLUMNS = {
            "Number_blocks":            "Number of blocks",
            "Average_coverage":         "Average synteny coverage (%)",
            "Coverage_min_genome_size": "Synteny coverage of smallest genome (%)",
            "N50_length":               "N50 length of synteny blocks (bp)",
        }

        with open(args.block_stats, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader)  # single summary row

        table_lines = [
            '<table id="tbl-blocks" class="data-table">',
            "<thead><tr><th>Statistic</th><th>Value</th></tr></thead>",
            "<tbody>",
        ]
        for col, label in BLOCK_STATS_COLUMNS.items():
            raw = row.get(col, "N/A")
            # Round to 2 decimal places if the value is a float
            try:
                value = f"{float(raw):,.2f}" if "." in raw else f"{int(raw):,}"
            except (ValueError, TypeError):
                value = raw
            table_lines.append(
                f"<tr><td>{html.escape(label)}</td><td>{html.escape(str(value))}</td></tr>"
            )
        table_lines.append("</tbody></table>")

        content = (
            f'<div class="table-wrap">{"".join(table_lines)}</div>'
            f'<p class="table-note">Source: {html.escape(args.block_stats)}</p>'
        )
    else:
        content = "<p><em>Synteny block stats not found.</em></p>"
    sections_html.append(section("Synteny block statistics", content, "synteny-blocks"))

# 3. Discontinuity reasons
    if args.discontinuity and os.path.exists(args.discontinuity):
        with open(args.discontinuity, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        table_lines = [
            '<table id="tbl-discontinuity" class="data-table">',
            "<thead><tr>"
            "<th>Synteny block discontinuity reason</th>"
            "<th>Count</th>"
            "</tr></thead>",
            "<tbody>",
        ]
        for row in rows:
            reason = html.escape(row["discontinuity_reason"])
            try:
                count = f"{int(row['discontinuity_reason_count']):,}"
            except (ValueError, TypeError):
                count = row["discontinuity_reason_count"]
            table_lines.append(
                f"<tr><td>{reason}</td><td>{count}</td></tr>"
            )
        table_lines.append("</tbody></table>")

        content = (
            f'<div class="table-wrap">{"".join(table_lines)}</div>'
            f'<p class="table-note">Source: {html.escape(args.discontinuity)}</p>'
        )
    else:
        content = "<p><em>Discontinuity reasons TSV not found.</em></p>"
    sections_html.append(section("Reasons for synteny block discontinuities", content, "discontinuity"))

    # 4. Mash divergence plot
    if args.mash_plot and os.path.exists(args.mash_plot):
        uri = encode_image(args.mash_plot)
        content = f"""
        <figure class="figure-wrap" style="max-width: 75%;">
          <img src="{uri}" alt="Mash divergence distributions">
          <figcaption>
            Pairwise mash distances across full-genome, syntenic, and
            non-syntenic regions. Wilcoxon signed-rank test p-value shown
            between syntenic and non-syntenic distributions.
            Image: {args.mash_plot}
          </figcaption>
        </figure>"""
    else:
        content = "<p><em>Mash divergence plot not found.</em></p>"
    sections_html.append(section("Mash divergence distributions", content, "divergence"))

# 5. Ribbon plot
    if args.ribbon_plot and os.path.exists(args.ribbon_plot):
        uri = encode_image(args.ribbon_plot)
        content = f"""
        <figure class="figure-wrap" style="max-width: 100%;">
          <img src="{uri}" alt="ntSynt-viz ribbon plot" style="width: 100%;">
          <figcaption>
            ntSynt-viz ribbon plot showing synteny blocks across all assemblies.
            Image: {args.ribbon_plot}
          </figcaption>
        </figure>"""
    else:
        content = "<p><em>Ribbon plot not found.</em></p>"
    sections_html.append(section("Synteny ribbon plot", content, "ribbon"))

    return HTML_TEMPLATE.format(
        group=html.escape(group_display),
        date=date_str,
        sections="\n".join(sections_html),
    )
PDF_CSS = """
@page {
    size: A4 portrait;
    margin: 1.5cm;
}

html {
    font-size: 11px;
}

main {
    max-width: 100%;
}

header, nav {
    page-break-after: avoid;
}

.report-section {
    page-break-inside: auto;
}

img {
    max-width: 100% !important;
}

/* Mash divergence plot — narrower since it's a portrait-style boxplot */
#divergence .figure-wrap {
    max-width: 100%;
}

#divergence .figure-wrap img {
    width: 100%;
}

/* Keep ribbon plot header and image on the same page */
#ribbon {
    page-break-inside: avoid;
}

/* Ribbon plot — full width since it's wide by nature */
#ribbon .figure-wrap {
    max-width: 100%;
    width: 100%;
    page-break-before: avoid;
    margin-left: -7.5%;
}

#ribbon .figure-wrap img {
    width: 100%;
}

"""


def save_pdf(html_content: str, pdf_path: str) -> None:
    """Render the HTML report to PDF using WeasyPrint."""
    try:
        from weasyprint import HTML, CSS
        from weasyprint.text.fonts import FontConfiguration
        font_config = FontConfiguration()
        html_obj = HTML(string=html_content)
        css = CSS(string=PDF_CSS, font_config=font_config)
        html_obj.write_pdf(pdf_path, stylesheets=[css], font_config=font_config)
    except ImportError:
        print(
            "WARNING: WeasyPrint not installed — skipping PDF output.\n"
            "Install with: mamba install -c conda-forge weasyprint",
            file=sys.stderr,
        )

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate a self-contained HTML report for an ntSynt run.")
    p.add_argument("--abyss-fac",     metavar="TSV",  help="abyss-fac summary TSV")
    p.add_argument("--block-stats",   metavar="TSV",  help="ntSynt synteny block stats TSV")
    p.add_argument("--ribbon-plot",   metavar="PNG",  help="ntSynt-viz ribbon plot PNG")
    p.add_argument("--discontinuity", metavar="TSV",  help="Block discontinuity reasons TSV")
    p.add_argument("--mash-plot",     metavar="PNG",  help="Mash divergence boxplot PNG")
    p.add_argument("--group",         required=True,  help="Taxonomic group name (used in title)")
    p.add_argument("--output",        required=True,  metavar="HTML", help="Output prefix file path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(args)
    html_path = Path(args.output).with_suffix(".html")
    pdf_path = str(Path(args.output).with_suffix(".pdf"))
    html_path.write_text(report, encoding="utf-8")
    save_pdf(report, pdf_path)
    print(f"Reports written to {str(html_path)} and {pdf_path}.")


if __name__ == "__main__":
    main()