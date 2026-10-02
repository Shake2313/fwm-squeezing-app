"""Build the static Site from its committed HTML source; no dependencies."""
from pathlib import Path
import argparse
import json
import re

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--refresh', action='store_true', help='Refresh from the parent standalone v3 atlas')
args = parser.parse_args()
source = ROOT / 'index.html'
if args.refresh:
    page = (ROOT.parent / 'fwm_quotient_structure_v3.html').read_text(encoding='utf-8')
    # Local source files are optional provenance, not deployment dependencies.
    # Preserve their names and all embedded references; link to the local catalog.
    page = re.sub(r'href="(?!#|https?:|data:)([^"\s]+)"', 'href="#sources"', page)
    page = page.replace('>원문 TeX</a>', '>원문 TeX 검색 키</a>')
    page = page.replace('>PDF</a>', '>원문 검색 키</a>')
    page = page.replace('>설계·변경 기록</a>', '>원문 근거</a>')
    page = page.replace('와 <a href="#sources">원문 검색 키</a>는 추가 검증용입니다.', '는 원문에서 대조할 때 사용합니다.')
    source.write_text(page, encoding='utf-8', newline='\n')

page = source.read_text(encoding='utf-8')
match = re.search(r'<script type="application/json" id="atlas-data">(.*?)</script>', page, re.S)
assert match, 'Missing embedded graph data'
graph = json.loads(match.group(1))
assert graph['version'] == 'v3'
assert sum(len(g['nodes']) for g in graph['groups']) == 61
assert len(graph['edges']) == 98
assert 'fwm_quotient_structure_v3.html' in page
assert not re.search(r'<(?:script|link)[^>]+(?:src|href)="https?://', page)
out = ROOT / 'out'
out.mkdir(exist_ok=True)
(out / 'index.html').write_text(page, encoding='utf-8', newline='\n')
print('Static Site built: 61 nodes, 98 relations, embedded offline explorer.')
