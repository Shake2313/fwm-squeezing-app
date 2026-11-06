"""Export current schematic/default inputs, without treating them as lab facts."""
from collections import Counter
from dataclasses import asdict
from pathlib import Path
import hashlib
import json
import subprocess
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT))
from sabes.layout import parts
from sabes.beamline import SetupSettings

def safe(value):
    if isinstance(value, float) and (value != value or abs(value) == float('inf')):
        return None
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [safe(v) for v in value]
    return value

cal = json.loads((ROOT / 'sabes/calibration_sim2025.json').read_text(encoding='utf-8'))
layouts = [v for v in vars(parts).values() if isinstance(v, parts.Layout)]
layouts = list({x.key: x for x in layouts}.values())
paths = ['sabes/layout/parts.py', 'sabes/layout/table.py', 'sabes/beamline.py', 'sabes/calibration.py', 'sabes/calibration_sim2025.json', 'sabes/detection.py']
hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
data = {
    'date': '2026-10-21', 'origin': 'code_snapshot', 'adoption_status': 'existing_implementation_only',
    'scope_note': 'AI export: display coordinates are arbitrary schematic units. Values and provenance labels are copied from code, not independently certified apparatus facts. Null defaults include NaN in source.',
    'git_branch': subprocess.check_output(['git', 'branch', '--show-current'], cwd=ROOT, text=True).strip(),
    'source_sha256': hashes,
    'setup_defaults': safe(asdict(SetupSettings())),
    'calibration': cal,
    'layouts': [asdict(x) for x in layouts],
    'routes': {k: list(v) for k, v in vars(parts).items() if k.startswith('ROUTE_')},
}
summary = {
    'date': data['date'], 'origin': data['origin'], 'adoption_status': data['adoption_status'],
    'coefficients': len(cal), 'coefficient_provenance': dict(Counter(x['provenance'] for x in cal.values())),
    'nodes': sum(len(x.nodes) for x in layouts), 'links': sum(len(x.links) for x in layouts),
    'layouts': [{ 'key': x.key, 'nodes': len(x.nodes), 'links': len(x.links), 'links_with_path_m': sum(y.path_m is not None for y in x.links)} for x in layouts],
    'source_sha256': hashes,
}
for name, content in [('sabes_implementation_inventory.json', data), ('sabes_inventory_summary.json', summary)]:
    (OUT / name).write_text(json.dumps(content, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False))
