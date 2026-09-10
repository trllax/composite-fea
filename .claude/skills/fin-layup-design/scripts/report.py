"""Fin layup candidate report -> a single self-contained HTML file.

  python report.py --baseline results/baseline --out report.html \
      --plybook-dir cases/fin_test_3/ax \
      results/ax/ax2r_uf0.66 results/ax/ax5r_uf0.60 ...

For each results dir it reads flex_kick.json + twist.json, embeds kick_shape.svg
(the deformed-blade side view, from run_tipweight.py --kick-bands N), draws the
ply stack from the matching <name>.csv in --plybook-dir, and prints the cut
sheet. The candidate table is sorted by flex. Deliberately plain HTML - this is
a work sheet, not an artifact; for a polished report hand the numbers to the
Artifact tool.
"""
from __future__ import annotations
import argparse, json, pathlib, html, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from cutsheet import cut_rows  # noqa: E402

MAT_COLOR = {  # material@angle -> fill
    ("twill_3k_198", 45): "#0f6f6c", ("twill_3k_198", 0): "#8fb6b4",
    ("hexcel_uni_231", 0): "#2b3033", ("hexcel_himax_biax_100", 45): "#5fa8c9",
    ("hexcel_himax_biax_100", 0): "#a9cede",
}


def _num(d, *path, default=None):
    for k in path:
        d = (d or {}).get(k)
    return d if d is not None else default


def load(run: pathlib.Path):
    fk = json.loads((run / "flex_kick.json").read_text())
    tw = {}
    if (run / "twist.json").exists():
        tw = json.loads((run / "twist.json").read_text())
    return dict(
        name=run.name, flex=_num(fk, "flex_n"),
        bucket=_num(fk, "kick_bucket"), frac=_num(fk, "kick_s_frac"),
        migr=_num(fk, "migration_delta_mm"),
        warn=_num(fk, "headline", "warning"),
        ktwist=_num(tw, "k_twist_nmm_per_rad"),
        theta=_num(tw, "theta_at_probe_deg"),
        twdeg=_num(tw, "phi_deg"), twwarn=_num(tw, "warning"),
        shape=(run / "kick_shape.svg"),
    )


def stack_svg(plybook: pathlib.Path) -> str:
    plies, _, _ = cut_rows(plybook)
    rows = list(csv_plies(plybook))
    w, x, gap = 460, 8, 2
    bar_h = 15
    out = [f'<svg viewBox="0 0 {w} {len(rows)*(bar_h+gap)+4}" '
           f'width="100%" style="max-width:520px">']
    for i, (zone, mat, ang, th) in enumerate(rows):
        y = 2 + i * (bar_h + gap)
        bw = (th / 0.30) * 300
        col = MAT_COLOR.get((mat, ang), "#999")
        out.append(f'<rect x="{x}" y="{y}" width="{bw:.0f}" height="{bar_h}" '
                   f'rx="2" fill="{col}"/>')
        out.append(f'<text x="{x+bw+8:.0f}" y="{y+bar_h-3}" '
                   f'font-family="ui-monospace,monospace" font-size="10.5" '
                   f'fill="currentColor">{i+1} · {zone} · '
                   f'{mat.split("_")[-1]}@{ang}° · {th:.3f}</text>')
    out.append("</svg>")
    return "\n".join(out)


def csv_plies(plybook: pathlib.Path):
    import csv
    txt = plybook.read_text().splitlines()
    for r in csv.DictReader(l for l in txt
                            if l.strip() and not l.lstrip().startswith("#")):
        yield (r["zone"], r["material"], int(float(r["angle_deg"])),
               float(r["thickness_mm"]))


def cutsheet_html(plybook: pathlib.Path) -> str:
    plies, per_mat, mass = cut_rows(plybook)
    r = ['<table class="cut"><thead><tr><th>fabric</th><th>gsm</th><th>plies</th>'
         '<th>orient</th><th>m²</th><th>+35% nest</th><th>roll</th></tr></thead><tbody>']
    for m, d in per_mat.items():
        from cutsheet import GSM
        ang = "/".join(f"{a}°" for a in sorted(d["angles"]))
        roll = (f"{d['roll_m']:.2f} m @ {d['roll_w_in']}\"" if d["roll_m"] else "—")
        r.append(f"<tr><td>{m}</td><td>{GSM[m]:.0f}</td><td>{d['plies']}</td>"
                 f"<td>{ang}</td><td>{d['fabric_m2']:.3f}</td>"
                 f"<td>{d['nest_m2']:.3f}</td><td>{roll}</td></tr>")
    r.append(f"</tbody></table><p class='mass'>dry fabric ~{mass:.0f} g → "
             f"~{mass*1.2:.0f} g cured</p>")
    return "".join(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", type=pathlib.Path)
    ap.add_argument("--baseline", type=pathlib.Path)
    ap.add_argument("--plybook-dir", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("fin_report.html"))
    ap.add_argument("--title", default="Fin layup candidates")
    a = ap.parse_args()

    cands = [load(r) for r in a.runs]
    base = load(a.baseline) if a.baseline else None
    kbase = base["ktwist"] if base and base["ktwist"] else None
    cands.sort(key=lambda c: (c["flex"] is None, c["flex"] or 0))

    def ratio(c):
        if not (kbase and c["ktwist"] and c["theta"] and base["theta"]):
            return "—"
        tag = "" if abs(c["theta"] - base["theta"]) <= 5 else " (θ off)"
        return f"{c['ktwist']/kbase:.2f}×{tag}"

    rowsrc = ([base] if base else []) + cands
    trs = []
    for c in rowsrc:
        cls = "base" if c is base else ""
        nm = c["name"].split("_uf")[0]
        trs.append(
            f"<tr class='{cls}'><td>{html.escape(nm)}</td>"
            f"<td>{c['flex']:.2f}</td><td>{c['bucket']}</td><td>{c['frac']:.2f}</td>"
            f"<td>{c['migr']:.1f}</td><td>{(c['ktwist'] or 0):.0f}</td>"
            f"<td>{(c['theta'] or 0):.0f}°</td>"
            f"<td>{'1.00×' if c is base else ratio(c)}</td>"
            f"<td>{c['warn'] or c['twwarn'] or ''}</td></tr>")

    blocks = []
    for c in cands:
        nm = c["name"].split("_uf")[0]
        pb = a.plybook_dir / f"{nm}.csv"
        shape = ""
        if c["shape"].exists():
            svg = c["shape"].read_text()
            svg = re.sub(r'width="[^"]*"', 'width="100%"', svg, count=1)
            shape = f"<figure class='shape'>{svg}<figcaption>deformed blade at "
            shape += "the 90° drive · colour = |curvature|</figcaption></figure>"
        stack = stack_svg(pb) if pb.exists() else "<p>(no ply book)</p>"
        cut = cutsheet_html(pb) if pb.exists() else ""
        blocks.append(f"""
        <section><h2>{html.escape(nm)}</h2>
        <p class='hd'>flex {c['flex']:.1f} N · kick {c['bucket']} f={c['frac']:.2f}
        · migration {c['migr']:.1f} mm · k_twist {(c['ktwist'] or 0):.0f}
        @ θ {(c['theta'] or 0):.0f}° · {ratio(c) if kbase else ''}</p>
        <div class='two'><div>{stack}</div><div>{shape}</div></div>
        {cut}</section>""")

    doc = f"""<!doctype html><meta charset=utf-8>
<title>{html.escape(a.title)}</title>
<style>
 :root{{color-scheme:light dark}}
 body{{font:15px/1.6 system-ui,sans-serif;max-width:1000px;margin:0 auto;padding:24px}}
 h1{{font-size:24px}} h2{{font-size:18px;margin-top:34px}}
 table{{border-collapse:collapse;width:100%;font:13px/1.5 ui-monospace,monospace;margin:14px 0}}
 th,td{{border-bottom:1px solid #8883;padding:6px 10px;text-align:right}}
 th:first-child,td:first-child{{text-align:left;font-family:system-ui}}
 tr.base td{{opacity:.7}}
 .hd{{font:13px ui-monospace,monospace;opacity:.8}}
 .two{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin:12px 0;align-items:start}}
 @media(max-width:700px){{.two{{grid-template-columns:1fr}}}}
 figure.shape{{margin:0}} figcaption{{font:11px ui-monospace,monospace;opacity:.7;margin-top:4px}}
 table.cut{{font-size:12px}} .mass{{font:12px ui-monospace,monospace;opacity:.8}}
</style>
<h1>{html.escape(a.title)}</h1>
<p class='hd'>generated by .claude/skills/fin-layup-design/scripts/report.py</p>
<table><thead><tr><th>design</th><th>flex N</th><th>kick</th><th>f</th>
<th>migr mm</th><th>k_twist</th><th>θ probe</th><th>× base</th><th>warn</th></tr></thead>
<tbody>{''.join(trs)}</tbody></table>
{''.join(blocks)}
"""
    a.out.write_text(doc, encoding="utf-8")
    print(f"wrote {a.out}  ({len(cands)} candidates)")


if __name__ == "__main__":
    main()
