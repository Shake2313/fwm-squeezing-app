"""Build a fully self-contained User's Guide for static serving.

Reads the editable guide in docs/, inlines every referenced PNG as a base64
data URI, and writes a single portable HTML file to static/ so the Streamlit app
can serve it at app/static/GABES_User_Guide.html — viewable from any computer,
no sibling asset files needed.

Run after editing docs/GABES_User_Guide_v2.html:
    python docs/build_static_guide.py
"""
import base64
import re
from pathlib import Path

DOCS = Path(__file__).resolve().parent
SRC = DOCS / "GABES_User_Guide_v2.html"
OUT = DOCS.parent / "static" / "GABES_User_Guide.html"


def _inline(match):
    rel = match.group(1)
    data = base64.b64encode((DOCS / rel).read_bytes()).decode("ascii")
    return f'src="data:image/png;base64,{data}"'


def main():
    html = SRC.read_text(encoding="utf-8")
    n_before = len(re.findall(r'src="userguide_assets/[^"]+\.png"', html))
    html = re.sub(r'src="(userguide_assets/[^"]+\.png)"', _inline, html)
    leftover = re.findall(r'src="(?!data:)[^"]+"', html)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"inlined {n_before} images -> {OUT}  ({kb:.0f} KB)")
    if leftover:
        print("WARNING: non-inlined external refs remain:", leftover)


if __name__ == "__main__":
    main()
