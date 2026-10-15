"""Summary statistics, control correction, figures and an HTML report."""
from __future__ import annotations

import base64
import html
import io
import math
from collections import Counter
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * max(0.0, c - h), 100 * min(1.0, c + h)


def metrics(c: Counter) -> Dict[str, float]:
    n = c.get("informative", 0)
    mod = n - c.get("unmodified", 0)
    lo, hi = wilson(mod, n)
    f = (lambda x: 100.0 * x / n) if n else (lambda x: float("nan"))
    return dict(informative_reads=n, modified_reads=mod, pct_modified=f(mod), pct_modified_ci95_lo=lo,
                pct_modified_ci95_hi=hi, pct_frameshift=f(c.get("frameshift", 0)),
                pct_inframe=f(c.get("in-frame", 0)), pct_exon_boundary=f(c.get("exon-boundary", 0)),
                pct_noncoding_indel=f(c.get("noncoding", 0)), pct_substitution_only=f(c.get("substitution", 0)))


def correct(t: float, c: float) -> float:
    """Background-correct an editing percentage with the unedited control:
    (T - C) / (100 - C) * 100, floored at 0 (Abbott-type correction: the control
    rate is the probability that an unedited molecule is *called* edited)."""
    if any(map(lambda x: x is None or (isinstance(x, float) and math.isnan(x)), (t, c))) or c >= 100:
        return float("nan")
    return max(0.0, (t - c) / (100 - c) * 100)


def _fig_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def plot_outcomes(rows: List[dict], path: str):
    names = [r["sample"] for r in rows]
    cats = [("pct_frameshift", "frameshift", "#b2182b"), ("pct_exon_boundary", "exon-boundary del", "#ef8a62"),
            ("pct_inframe", "in-frame", "#67a9cf"), ("pct_noncoding_indel", "non-coding indel", "#bbbbbb")]
    fig, ax = plt.subplots(figsize=(max(4, 1.1 * len(names) + 2), 3.6))
    bottom = [0.0] * len(names)
    for key, lab, col in cats:
        vals = [r.get(key, 0) or 0 for r in rows]
        ax.bar(names, vals, bottom=bottom, label=lab, color=col, edgecolor="white", linewidth=0.5)
        bottom = [b + v for b, v in zip(bottom, vals)]
    for i, r in enumerate(rows):
        if r.get("pct_modified_ci95_hi") == r.get("pct_modified_ci95_hi"):
            ax.errorbar(i, r["pct_modified"], yerr=[[r["pct_modified"] - r["pct_modified_ci95_lo"]],
                                                     [r["pct_modified_ci95_hi"] - r["pct_modified"]]],
                        color="k", capsize=3, lw=1)
    ax.set_ylabel("% of informative reads")
    ax.set_ylim(0, 100)
    ax.set_title("Indel outcomes at the cut site")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def plot_sizes(sample: str, allele_rows: List[tuple], path: str, lo: int = -40, hi: int = 15):
    c = Counter()
    tot = sum(r[3] for r in allele_rows)
    for lab, cat, net, n, _ in allele_rows:
        if cat == "unmodified":
            continue
        c[max(lo, min(hi, net))] += n
    xs = list(range(lo, hi + 1))
    ys = [100 * c[x] / tot if tot else 0 for x in xs]
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    cols = ["#b2182b" if x % 3 else "#2166ac" for x in xs]
    ax.bar(xs, ys, color=cols, width=0.85)
    ax.set_xlabel(f"net indel size (bp); bins at {lo} / +{hi} include larger events")
    ax.set_ylabel("% of reads")
    ax.set_title(f"{sample}: indel size distribution (red = frameshift, blue = in-frame)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def html_report(path: str, title: str, preamble: List[str], summary_rows: List[dict], columns: List[str],
                figs: List[tuple], allele_tables: Dict[str, List[tuple]], notes: List[str]):
    def table(rows, cols):
        h = "<table><tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in cols) + "</tr>"
        for r in rows:
            cells = []
            for c in cols:
                v = r.get(c, "") if isinstance(r, dict) else r[cols.index(c)]
                if isinstance(v, float):
                    v = "NA" if math.isnan(v) else f"{v:.2f}"
                cells.append(f"<td>{html.escape(str(v))}</td>")
            h += "<tr>" + "".join(cells) + "</tr>"
        return h + "</table>"
    parts = [f"<html><head><meta charset='utf-8'><title>{html.escape(title)}</title><style>"
             "body{font-family:Helvetica,Arial,sans-serif;max-width:1100px;margin:2em auto;color:#222}"
             "table{border-collapse:collapse;font-size:12px;margin:0.5em 0 1.5em}"
             "td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}th{background:#f2f2f2}"
             "pre{background:#f7f7f7;padding:8px;font-size:11px;overflow-x:auto}"
             "code{font-size:11px}</style></head><body>", f"<h1>{html.escape(title)}</h1>"]
    parts.append("<pre>" + html.escape("\n".join(preamble)) + "</pre>")
    parts.append("<h2>Summary</h2>" + table(summary_rows, columns))
    for cap, fig in figs:
        parts.append(f"<h3>{html.escape(cap)}</h3><img style='max-width:100%' src='data:image/png;base64,{_fig_b64(fig)}'>")
    for s, rows in allele_tables.items():
        parts.append(f"<h3>{html.escape(s)}: top alleles (window used for the headline call)</h3>")
        parts.append(table([dict(allele=a, category=c, net=n, reads=k, pct=p) for a, c, n, k, p in rows],
                           ["allele", "category", "net", "reads", "pct"]))
    parts.append("<h2>Notes</h2><ul>" + "".join(f"<li>{html.escape(n)}</li>" for n in notes) + "</ul>")
    parts.append("</body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(parts))
