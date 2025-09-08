#!/usr/bin/env python3
"""
Interactive dotplot from a nucmer .delta file using plotly.
# requires mummer, plotly, and python3
- --ref_contigs and --qry_contigs let you include only specific contigs

Usage:
  python dotplot_from_delta.py alignment.delta \
      --min_len 1000 --min_idy 90 \
      --ref_contigs chr2L chr2R \
      --qry_contigs scaffold123
"""
import argparse
import os
import shutil
import subprocess
import sys

import plotly.graph_objects as go
from plotly.offline import plot as plot_offline

PLOT_POINTS = False  # set True to plot midpoints as dots instead of segments


def run_show_coords(delta_path):
    exe = shutil.which("show-coords")
    if not exe:
        sys.stderr.write(
            "ERROR: `show-coords` not found. Please install mummer\n"
        )
        sys.exit(1)
    # -r: sort by reference, -c: include percent coverage, -l: long format (names),
    # -T: tab-delimited, -H: no header
    cmd = [exe, "-rclTH", delta_path]
    try:
        out = subprocess.check_output(cmd, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"ERROR running {' '.join(cmd)}:\n{e}\n")
        sys.exit(1)
    return out.strip().splitlines()


def parse_coords(lines):
    """
    extract ref and qry names from the remaining fields
    """
    rows = []
    for ln in lines:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 7:
            parts = ln.split()
        try:
            S1 = int(parts[0]); E1 = int(parts[1]); S2 = int(parts[2]); E2 = int(parts[3])
            LEN1 = int(parts[4]); LEN2 = int(parts[5]); IDY = float(parts[6])
        except (ValueError, IndexError):
            continue

        ref_name = None
        qry_name = None
        tail = [p for p in parts[7:] if any(c.isalpha() for c in p)]
        if len(tail) >= 2:
            ref_name = tail[-2]
            qry_name = tail[-1]

        rows.append({
            "S1": S1, "E1": E1, "S2": S2, "E2": E2,
            "LEN1": LEN1, "LEN2": LEN2, "IDY": IDY,
            "ref": ref_name or "ref", "qry": qry_name or "qry"
        })
    return rows


def build_traces(rows, min_len, min_idy, ref_filter=None, qry_filter=None):
    """
    ref_filter, qry_filter: optional sets of contig names to include
    """
    fwd_x = []; fwd_y = []; fwd_text = []
    rev_x = []; rev_y = []; rev_text = []
    refs = set(); qrys = set()

    for r in rows:
        if ref_filter and r["ref"] not in ref_filter:
            continue
        if qry_filter and r["qry"] not in qry_filter:
            continue
        if r["LEN1"] < min_len and r["LEN2"] < min_len:
            continue
        if r["IDY"] < min_idy:
            continue

        refs.add(r["ref"]); qrys.add(r["qry"])
        ori = (r["E1"] - r["S1"]) * (r["E2"] - r["S2"])

        x0, x1 = r["S1"], r["E1"]
        y0, y1 = r["S2"], r["E2"]

        hover = (
            f"<b>{r['ref']}</b>: {r['S1']:,}-{r['E1']:,} (len1={r['LEN1']:,})<br>"
            f"<b>{r['qry']}</b>: {r['S2']:,}-{r['E2']:,} (len2={r['LEN2']:,})<br>"
            f"IDY: {r['IDY']:.2f}%"
        )

        if PLOT_POINTS:
            xm = (x0 + x1) / 2.0
            ym = (y0 + y1) / 2.0
            if ori >= 0:
                fwd_x.append(xm); fwd_y.append(ym); fwd_text.append(hover)
            else:
                rev_x.append(xm); rev_y.append(ym); rev_text.append(hover)
        else:
            if ori >= 0:
                fwd_x += [x0, x1, None]; fwd_y += [y0, y1, None]; fwd_text += [hover, hover, hover]
            else:
                rev_x += [x0, x1, None]; rev_y += [y0, y1, None]; rev_text += [hover, hover, hover]

    fwd = go.Scattergl(
        x=fwd_x, y=fwd_y, mode="lines" if not PLOT_POINTS else "markers",
        name="Forward", hoverinfo="text", text=fwd_text, line=dict(width=1)
    )
    rev = go.Scattergl(
        x=rev_x, y=rev_y, mode="lines" if not PLOT_POINTS else "markers",
        name="Reverse", hoverinfo="text", text=rev_text, line=dict(width=1, dash="dot")
    )
    return fwd, rev, sorted(refs), sorted(qrys)


def build_layout(ref_label="Reference", qry_label="Query"):
    return go.Layout(
        title=f"Interactive Dotplot ({ref_label} vs {qry_label})",
        xaxis=dict(title=f"{ref_label} position", zeroline=False, rangemode="tozero"),
        yaxis=dict(title=f"{qry_label} position", zeroline=False, scaleanchor=None),
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=60, r=10, t=60, b=60),
    )


def make_filters(fig, refs, qrys):
    if len(refs) <= 1 and len(qrys) <= 1:
        return fig
    fig.update_layout(
        annotations=[
            dict(
                text=f"Ref contigs: {', '.join(refs[:20])}{'…' if len(refs)>20 else ''}<br>"
                     f"Qry contigs: {', '.join(qrys[:20])}{'…' if len(qrys)>20 else ''}",
                xref="paper", yref="paper", x=0, y=1.12, showarrow=False, align="left",
                font=dict(size=10)
            )
        ]
    )
    return fig


def main():
    ap = argparse.ArgumentParser(description="Interactive dotplot from a nucmer .delta file")
    ap.add_argument("delta", help="Path to nucmer .delta")
    ap.add_argument("--min_len", type=int, default=1000, help="Min alignment length to plot (default 1000)")
    ap.add_argument("--min_idy", type=float, default=0.0, help="Min percent identity to plot (default 0)")
    ap.add_argument("--out", default=None, help="Output HTML file (default: <delta>.dotplot.html)")
    ap.add_argument("--ref_contigs", nargs="+", help="Reference contigs to include (space-separated)")
    ap.add_argument("--qry_contigs", nargs="+", help="Query contigs to include (space-separated)")
    args = ap.parse_args()

    if not os.path.isfile(args.delta):
        sys.stderr.write(f"ERROR: File not found: {args.delta}\n")
        sys.exit(1)

    lines = run_show_coords(args.delta)
    rows = parse_coords(lines)
    if not rows:
        sys.stderr.write("No alignments parsed. Check your .delta and show-coords version/options.\n")
        sys.exit(1)

    ref_filter = set(args.ref_contigs) if args.ref_contigs else None
    qry_filter = set(args.qry_contigs) if args.qry_contigs else None

    fwd, rev, refs, qrys = build_traces(
        rows, args.min_len, args.min_idy, ref_filter=ref_filter, qry_filter=qry_filter
    )
    layout = build_layout("Reference", "Query")
    fig = go.Figure(data=[fwd, rev], layout=layout)
    fig = make_filters(fig, refs, qrys)

    out_html = args.out or os.path.splitext(args.delta)[0] + ".dotplot.html"
    plot_offline(fig, filename=out_html, auto_open=True)
    print(f"Wrote {out_html}")


if __name__ == "__main__":
    main()
