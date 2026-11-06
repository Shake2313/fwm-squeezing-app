"""Read-only text inventory of the apparatus source archive; no Office COM."""
from pathlib import Path
from html.parser import HTMLParser
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parent
LAB = Path('G:/내 드라이브/원자광학연구실/Lab meeting')
SIM = Path('G:/내 드라이브/Docs/Gisung Sim')
NS = {'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}

class PlainHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.text = []
        self.hidden = 0
    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'):
            self.hidden += 1
    def handle_endtag(self, tag):
        if tag in ('script', 'style'):
            self.hidden = max(0, self.hidden - 1)
    def handle_data(self, data):
        if not self.hidden and data.strip():
            self.text.append(data.strip())

def scan(path):
    raw = path.read_bytes()
    item = {'path': str(path), 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(), 'sections': []}
    if path.suffix.lower() == '.pptx':
        with zipfile.ZipFile(path) as z:
            for stem in ('ppt/slides/slide', 'ppt/notesSlides/notesSlide'):
                names = [n for n in z.namelist() if re.fullmatch(re.escape(stem) + r'\d+\.xml', n)]
                for name in sorted(names, key=lambda n: int(re.search(r'(\d+)\.xml$', n).group(1))):
                    root = ET.fromstring(z.read(name))
                    item['sections'].append({'part': name, 'text': '\n'.join(t.text or '' for t in root.findall('.//a:t', NS))})
    elif path.suffix.lower() == '.html':
        parser = PlainHTML()
        parser.feed(raw.decode('utf-8-sig'))
        item['sections'].append({'part': 'html text', 'text': '\n'.join(parser.text)})
    else:
        item['sections'].append({'part': 'text', 'text': raw.decode('utf-8-sig')})
    return item

paths = list(LAB.glob('*.pptx')) + list(LAB.glob('*.html'))
paths += list((LAB / 'lab-context').glob('*.md'))
paths += list(LAB.glob('2026*.md'))
paths += list(SIM.glob('*심기성*.pptx'))
paths += [LAB / 'Media & Reference/20260521/Etalon.md']
paths += [SIM.parent / 'ECDL laser 구성요소 및 제작방법_2024.08.02.pptx']
paths += list((LAB / 'Media & Reference/Templetes_PPTX').glob('*Squeezing*'))
records = []
for path in sorted(set(paths)):
    if not path.is_file() or path.name.startswith('~$'):
        continue
    try:
        records.append(scan(path))
    except (zipfile.BadZipFile, UnicodeError, ET.ParseError) as exc:
        records.append({'path': str(path), 'error': str(exc), 'sections': []})
(OUT / 'lab_sources_text.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
manifest = [{k: v for k, v in r.items() if k != 'sections'} | {'section_count': len(r['sections'])} for r in records]
(OUT / 'lab_sources_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'files': len(records), 'pptx': sum(r['path'].endswith('.pptx') for r in records), 'sections': sum(len(r['sections']) for r in records)}, ensure_ascii=False))
