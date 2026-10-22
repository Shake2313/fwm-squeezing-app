"""Run the required regression and preserve its completion/source receipt."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def source_state():
    files = [ROOT / "sabes_page.py", ROOT / "streamlit_app.py"]
    for name in ("gabes", "gabes_ui", "sabes", "tests"):
        files.extend((ROOT / name).rglob("*.py"))
        if name == "sabes":
            files.extend((ROOT / name).glob("*.json"))
    return {str(path.relative_to(ROOT)).replace("\\", "/"):
            hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(set(files))}


if __name__ == "__main__":
    before = source_state()
    (OUT / "validation_source_before.json").write_text(
        json.dumps(before, indent=2), encoding="utf-8")
    started = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    with (OUT / "full_pytest.log").open("w", encoding="utf-8") as log:
        result = subprocess.run([sys.executable, "-m", "pytest", "-q"],
                                cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    after = source_state()
    changed = [key for key in sorted(set(before) | set(after))
               if before.get(key) != after.get(key)]
    receipt = {"command": "python -m pytest -q", "started_at_utc": started,
               "finished_at_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": round(time.monotonic() - start, 2),
               "exit_code": result.returncode, "source_count": len(before),
               "changed_sources": changed,
               "log": "full_pytest.log"}
    (OUT / "full_pytest_receipt.json").write_text(
        json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps(receipt), flush=True)
    sys.exit(result.returncode)
