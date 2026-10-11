"""Turn paper/REPORT.md into paper/REPORT.html and paper/REPORT.pdf (headless Edge). A small converter, enough for this report."""
import html
import re
import subprocess
from pathlib import Path

PAPER = Path("D:/AAN/paper")
EDGE = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"

CSS = """
body{font-family:Georgia,'Times New Roman',serif;max-width:760px;margin:30px auto;padding:0 18px;color:#1c1c1c;line-height:1.5;font-size:11pt}
h1{font-size:20pt;line-height:1.2;margin-bottom:4px} h2{font-size:15pt;margin-top:26px;border-bottom:1px solid #ccc;padding-bottom:3px}
h3{font-size:12.5pt;margin-top:18px} p{margin:8px 0;text-align:left}
img{max-width:100%;display:block;margin:10px auto 2px} em.cap{display:block;font-size:9.5pt;color:#444;margin:2px 0 14px}
table{border-collapse:collapse;font-size:8.5pt;margin:8px 0 14px;width:100%} th,td{border:1px solid #bbb;padding:3px 5px;text-align:left} th{background:#eef2f5}
ul,ol{margin:6px 0 6px 20px} li{margin:3px 0} hr{border:0;border-top:1px solid #ddd;margin:20px 0}
@page{size:Letter;margin:0.7in}
"""


def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*(?!\s)(.+?)(?<!\s)\*(?!\*)", r"<em>\1</em>", s)
    s = re.sub(r"\^(-?\d+)", r"<sup>\1</sup>", s)
    return s


def convert(md):
    out, lines, i = [], md.split("\n"), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("!["):
            m = re.match(r"!\[(.*?)\]\((.*?)\)", ln)
            alt, src = m.group(1), m.group(2)
            out.append(f'<img src="{src}" alt="{html.escape(alt)}">')
            if re.match(r"Figure \d+\.", alt):
                out.append(f'<em class="cap">{inline(alt)}</em>')
        elif ln.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(lines[i]); i += 1
            i -= 1
            cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows if not re.match(r"^\|[-| ]+\|$", r)]
            t = "<table><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells[0]) + "</tr>"
            for r in cells[1:]:
                t += "<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>"
            out.append(t + "</table>")
        elif ln.startswith("### "):
            out.append(f"<h3>{inline(ln[4:])}</h3>")
        elif ln.startswith("## "):
            out.append(f"<h2>{inline(ln[3:])}</h2>")
        elif ln.startswith("# "):
            out.append(f"<h1>{inline(ln[2:])}</h1>")
        elif ln.strip() == "---":
            out.append("<hr>")
        elif re.match(r"^(- |\d+\. )", ln):
            ordered = bool(re.match(r"^\d+\. ", ln))
            tag = "ol" if ordered else "ul"
            items = []
            while i < len(lines) and re.match(r"^(- |\d+\. )", lines[i]):
                items.append(re.sub(r"^(- |\d+\. )", "", lines[i])); i += 1
            i -= 1
            out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>")
        elif ln.strip():
            out.append(f"<p>{inline(ln)}</p>")
        i += 1
    return "\n".join(out)


def main():
    md = (PAPER / "REPORT.md").read_text(encoding="utf-8")
    body = convert(md)
    doc = f"<!doctype html><html><head><meta charset='utf-8'><title>ME/CFS cell-type report</title><style>{CSS}</style></head><body>{body}</body></html>"
    (PAPER / "REPORT.html").write_text(doc, encoding="utf-8")
    pdf = PAPER / "REPORT.pdf"
    subprocess.run([EDGE, "--headless", "--disable-gpu", "--user-data-dir=D:/AAN_data/edge_tmp", "--no-pdf-header-footer", f"--print-to-pdf={pdf}",
                    (PAPER / "REPORT.html").as_uri()], check=True, timeout=180)
    print("wrote", pdf, pdf.stat().st_size // 1024, "KB")


if __name__ == "__main__":
    main()
