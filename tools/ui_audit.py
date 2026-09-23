"""
GABES UI audit — layout metrics and first-screen screenshots per scheme.

Drives a headless Chromium-family browser (Edge / Chrome) over the DevTools
protocol with tornado's websocket client, which Streamlit already depends on,
so the audit needs no extra packages. The UI redesign (docs/ui_redesign/)
records one run per phase and compares them.

Usage
    python tools/ui_audit.py --tag baseline            # spawns Streamlit on a free port
    python tools/ui_audit.py --tag p1 --url http://localhost:8501
    python tools/ui_audit.py --compare baseline p1     # markdown diff table

DOM contract the redesign must keep (or update here in the same commit):
  * stApp carries Streamlit's data-test-script-state attribute;
  * the scheme switcher is a role=combobox whose aria-label contains "Scheme",
    reachable with the sidebar closed (top bar), and its options are
    role=option elements whose text is the scheme title;
  * the main plot is the largest img / iframe / canvas inside stMain;
  * the control rail is Streamlit's sidebar (stSidebar / stSidebarContent);
  * the advanced toggle is a sidebar summary/button whose text contains
    "Advanced controls", "Show advanced" or "Hide advanced".
"""
from __future__ import annotations

import argparse
import asyncio
import base64
import itertools
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUDIT_ROOT = REPO / "docs" / "ui_redesign" / "audit"
DESKTOP = dict(width=1440, height=900, deviceScaleFactor=1, mobile=False)
MOBILE = dict(width=390, height=844, deviceScaleFactor=1, mobile=True)

BROWSER_CANDIDATES = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
]
BROWSER_NAMES = ["msedge", "google-chrome", "chromium", "chromium-browser", "chrome"]


def log(msg):
    print(f"[ui_audit] {msg}", file=sys.stderr, flush=True)


# ----------------------------------------------------------------------
# Scheme inventory (DOM-independent denominators)
# ----------------------------------------------------------------------
def _visible_by_default(spec, defaults):
    if getattr(spec, "hidden", False):
        return False
    for pname, allowed in (getattr(spec, "visible_if", None) or {}).items():
        cur = defaults.get(pname)
        if isinstance(allowed, (set, tuple, list)):
            if cur not in allowed:
                return False
        elif cur != allowed:
            return False
    return True


def scheme_inventory():
    sys.path.insert(0, str(REPO))
    from gabes import schemes

    rows = []
    for scheme in schemes.all_schemes():
        specs = scheme.param_schema()
        defaults = {sp.name: sp.default for sp in specs}
        visible = [sp for sp in specs if _visible_by_default(sp, defaults)]
        rows.append(dict(
            name=scheme.name,
            title=scheme.title,
            n_main_controls=sum(1 for sp in visible if not sp.advanced),
            n_advanced_controls=sum(1 for sp in visible if sp.advanced),
        ))
    return rows


# ----------------------------------------------------------------------
# Streamlit server + browser processes
# ----------------------------------------------------------------------
def _free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _http_ok(url, timeout=2.0):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.status == 200
    except Exception:
        return False


def spawn_streamlit():
    port = _free_port()
    cmd = [sys.executable, "-m", "streamlit", "run", "streamlit_app.py",
           "--server.port", str(port), "--server.headless", "true",
           "--browser.gatherUsageStats", "false"]
    proc = subprocess.Popen(cmd, cwd=str(REPO), stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    url = f"http://127.0.0.1:{port}"
    deadline = time.time() + 120
    while time.time() < deadline:
        if proc.poll() is not None:
            raise RuntimeError("streamlit exited during startup")
        if _http_ok(url + "/_stcore/health"):
            log(f"streamlit up at {url}")
            return proc, url
        time.sleep(0.5)
    proc.terminate()
    raise RuntimeError("streamlit did not become healthy within 120 s")


def find_browser(explicit=None):
    for cand in [explicit, os.environ.get("GABES_AUDIT_BROWSER"), *BROWSER_CANDIDATES]:
        if cand and Path(cand).exists():
            return cand
    for name in BROWSER_NAMES:
        found = shutil.which(name)
        if found:
            return found
    raise RuntimeError("no Edge/Chrome found; set GABES_AUDIT_BROWSER")


def launch_browser(executable):
    profile = tempfile.mkdtemp(prefix="gabes_audit_")
    args = [executable, "--headless=new", "--disable-gpu", "--no-first-run",
            "--no-default-browser-check", "--hide-scrollbars",
            "--remote-debugging-port=0", "--remote-allow-origins=*",
            f"--user-data-dir={profile}", "--window-size=1440,900", "about:blank"]
    proc = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    port_file = Path(profile) / "DevToolsActivePort"
    deadline = time.time() + 30
    while time.time() < deadline:
        if port_file.exists():
            text = port_file.read_text(encoding="utf-8").split()
            if text:
                return proc, profile, int(text[0])
        time.sleep(0.2)
    proc.terminate()
    raise RuntimeError("browser did not expose a DevTools port")


def page_ws_url(port):
    with urllib.request.urlopen(f"http://127.0.0.1:{port}/json/list", timeout=5) as resp:
        targets = json.loads(resp.read().decode("utf-8"))
    for target in targets:
        if target.get("type") == "page":
            return target["webSocketDebuggerUrl"]
    req = urllib.request.Request(f"http://127.0.0.1:{port}/json/new?about:blank", method="PUT")
    with urllib.request.urlopen(req, timeout=5) as resp:
        return json.loads(resp.read().decode("utf-8"))["webSocketDebuggerUrl"]


# ----------------------------------------------------------------------
# Minimal DevTools client
# ----------------------------------------------------------------------
class CDP:
    def __init__(self, conn):
        self.conn = conn
        self.ids = itertools.count(1)
        self.pending = {}
        self.reader = asyncio.ensure_future(self._read())

    @classmethod
    async def connect(cls, ws_url):
        from tornado.websocket import websocket_connect
        conn = await websocket_connect(ws_url, max_message_size=64 * 1024 * 1024)
        return cls(conn)

    async def _read(self):
        while True:
            msg = await self.conn.read_message()
            if msg is None:
                for fut in self.pending.values():
                    if not fut.done():
                        fut.set_exception(RuntimeError("DevTools connection closed"))
                return
            data = json.loads(msg)
            fut = self.pending.pop(data.get("id"), None)
            if fut is not None and not fut.done():
                if "error" in data:
                    fut.set_exception(RuntimeError(str(data["error"])))
                else:
                    fut.set_result(data.get("result", {}))

    async def send(self, method, **params):
        msg_id = next(self.ids)
        fut = asyncio.get_running_loop().create_future()
        self.pending[msg_id] = fut
        await self.conn.write_message(json.dumps({"id": msg_id, "method": method, "params": params}))
        return await asyncio.wait_for(fut, timeout=120)

    async def js(self, expression):
        res = await self.send("Runtime.evaluate", expression=expression,
                              returnByValue=True, awaitPromise=True)
        if "exceptionDetails" in res:
            raise RuntimeError(res["exceptionDetails"].get("text", "JS error"))
        return res.get("result", {}).get("value")

    async def click(self, x, y):
        for kind in ("mouseMoved", "mousePressed", "mouseReleased"):
            await self.send("Input.dispatchMouseEvent", type=kind, x=x, y=y,
                            button="left", clickCount=1)

    async def close(self):
        self.reader.cancel()
        self.conn.close()


# ----------------------------------------------------------------------
# Page-side probes
# ----------------------------------------------------------------------
JS_READY = r"""
(() => {
  const app = document.querySelector('[data-testid="stApp"]');
  if (!app) return {ready: false, state: 'noapp'};
  const state = app.getAttribute('data-test-script-state');
  const main = document.querySelector('[data-testid="stMain"]');
  const big = main ? [...main.querySelectorAll('img,iframe,canvas')].filter(e => {
    const r = e.getBoundingClientRect(); return r.width * r.height > 40000; }) : [];
  const spinning = !!document.querySelector('[data-testid="stSpinner"]');
  const sw = [...document.querySelectorAll('[role="combobox"]')]
    .find(e => (e.getAttribute('aria-label') || '').includes('Scheme'));
  return {ready: state === 'notRunning' && !spinning && big.length > 0,
          state, plots: big.length, scheme: sw ? sw.getAttribute('aria-label') : null,
          sig: big.map(e => { const r = e.getBoundingClientRect();
                              return [Math.round(r.top), Math.round(r.height)].join(':'); }).join('|')};
})()
"""

JS_MEASURE = r"""
(() => {
  const vh = innerHeight, vw = innerWidth;
  const main = document.querySelector('[data-testid="stMain"]');
  if (main) main.scrollTop = 0;
  const rail = document.querySelector('[data-testid="stSidebarContent"]')
            || document.querySelector('[data-testid="stSidebar"]');
  const sidebar = document.querySelector('[data-testid="stSidebar"]');
  const cands = main ? [...main.querySelectorAll('img,iframe,canvas')]
    .map(e => ({e, r: e.getBoundingClientRect()}))
    .filter(o => o.r.width > 0 && o.r.height > 0) : [];
  cands.sort((a, b) => b.r.width * b.r.height - a.r.width * a.r.height);
  let plot = null;
  if (cands.length) {
    const top = cands[0];
    plot = {top: top.r.top, bottom: top.r.bottom, width: top.r.width, height: top.r.height,
            kind: top.e.tagName.toLowerCase()};
    if (top.e.tagName === 'IFRAME') {
      try {
        const imgs = [...top.e.contentDocument.querySelectorAll('img')]
          .map(i => i.getBoundingClientRect()).filter(r => r.width > 0)
          .sort((a, b) => b.width * b.height - a.width * a.height);
        if (imgs.length) plot = {top: top.r.top + imgs[0].top, bottom: top.r.top + imgs[0].bottom,
                                 width: imgs[0].width, height: imgs[0].height, kind: 'iframe>img'};
      } catch (_e) {}
    }
    for (const k of ['top', 'bottom', 'width', 'height']) plot[k] = Math.round(plot[k]);
    plot.fully_visible = plot.top >= 0 && plot.bottom <= vh + 1;
  }
  const leaves = [...document.querySelectorAll('[data-testid="stMain"] *, [data-testid="stSidebar"] *')]
    .filter(e => {
      if (e.children.length) return false;
      if (!(e.textContent || '').trim()) return false;
      const r = e.getBoundingClientRect();
      return r.width > 0 && r.height > 0 && r.top < vh && r.bottom > 0 && r.left < vw && r.right > 0;
    }).length;
  const railOpen = !!sidebar && sidebar.getAttribute('aria-expanded') !== 'false'
                   && sidebar.getBoundingClientRect().width > 60;
  return {
    viewport: [vw, vh],
    plot,
    main_scroll_height: main ? main.scrollHeight : null,
    rail_open: railOpen,
    rail_scroll_height: rail && railOpen ? rail.scrollHeight : null,
    // scrollHeight never drops below the viewport, so a rail that now fits on
    // one screen reads as a flat 900 px. This is what the controls occupy.
    rail_content_height: (() => {
      const box = rail && railOpen
        ? rail.querySelector('[data-testid="stSidebarUserContent"]') : null;
      if (!box) return null;
      const kids = [...box.querySelectorAll(
        ':scope > div > [data-testid="stVerticalBlock"] > *')];
      if (!kids.length) return Math.round(box.getBoundingClientRect().height);
      const top = kids[0].getBoundingClientRect().top;
      const bottom = kids[kids.length - 1].getBoundingClientRect().bottom;
      return Math.round(bottom - top);
    })(),
    help_glyphs: document.querySelectorAll(
      '[data-testid="stTooltipIcon"],[data-testid="stTooltipHoverTarget"]').length,
    // The "?" next to a widget label. The same testid also wraps buttons that
    // merely carry a tooltip (About, SABES), which draw no glyph at all.
    help_icons: [...document.querySelectorAll('[data-testid="stTooltipIcon"]')]
      .filter(e => !e.querySelector('button')).length,
    iframes: document.querySelectorAll('iframe').length,
    main_expanders: main ? main.querySelectorAll('[data-testid="stExpander"]').length : 0,
    text_leaves_first_viewport: leaves,
    dom_nodes: document.querySelectorAll('*').length,
  };
})()
"""

JS_TRANSFER = r"""
(() => {
  const nav = performance.getEntriesByType('navigation')[0];
  const res = performance.getEntriesByType('resource');
  const bytes = (nav ? nav.transferSize : 0) + res.reduce((s, e) => s + (e.transferSize || 0), 0);
  return {transfer_kb: Math.round(bytes / 1024), requests: res.length + 1};
})()
"""

JS_ADVANCED_TOGGLE = r"""
(() => {
  const sb = document.querySelector('[data-testid="stSidebar"]');
  if (!sb) return false;
  // Streamlit's expander summary text starts with an icon ligature
  // ("keyboard_arrow_right"), so match anywhere and keep the tightest element.
  const hits = [...sb.querySelectorAll('summary,button')]
    .filter(e => /Advanced controls|Show advanced|Hide advanced/i.test(e.textContent || ''))
    .sort((a, b) => (a.textContent || '').length - (b.textContent || '').length);
  const el = hits[0];
  if (!el) return false;
  el.click();
  return true;
})()
"""

JS_SIDEBAR = r"""
((want) => {
  const sb = document.querySelector('[data-testid="stSidebar"]');
  const open = !!sb && sb.getAttribute('aria-expanded') !== 'false' && sb.getBoundingClientRect().width > 60;
  if (open === want) return 'ok';
  const sel = ['[data-testid="stSidebarCollapseButton"] button', '[data-testid="stExpandSidebarButton"]',
               '[data-testid="stSidebarCollapsedControl"] button',
               '[data-testid="stSidebarHeader"] button', '[data-testid="stSidebar"] [data-testid="stBaseButton-headerNoPadding"]'];
  for (const s of sel) { const b = document.querySelector(s); if (b) { b.click(); return 'clicked'; } }
  return 'missing';
})
"""


def _center_js(selector_expr):
    return (
        "(() => { const el = " + selector_expr + "; if (!el) return null;"
        " el.scrollIntoView({block: 'center'}); const r = el.getBoundingClientRect();"
        " return [r.left + r.width / 2, r.top + r.height / 2]; })()"
    )


async def wait_ready(cdp, *, expect_title=None, timeout=600.0, min_wait=0.0):
    start = time.time()
    last_sig = None
    stable = 0
    await asyncio.sleep(min_wait)
    while time.time() - start < timeout:
        try:
            st = await cdp.js(JS_READY)
        except Exception:
            st = {"ready": False}
        title_ok = expect_title is None or (st.get("scheme") or "").find(expect_title) >= 0
        if st.get("ready") and title_ok:
            stable = stable + 1 if st.get("sig") == last_sig else 1
            last_sig = st.get("sig")
            if stable >= 2:
                return round((time.time() - start) * 1000)
        else:
            stable = 0
        await asyncio.sleep(0.7)
    raise TimeoutError(f"page not ready after {timeout:.0f} s (expect {expect_title!r})")


async def switch_scheme(cdp, title):
    """Open the scheme switcher and pick `title`; retried because a rerun that is
    still settling can re-render the select and close its menu."""
    title_js = json.dumps(title)
    for _attempt in range(4):
        await wait_ready(cdp, min_wait=0.5)
        pos = await cdp.js(_center_js(
            "[...document.querySelectorAll('[role=\"combobox\"]')]"
            ".find(e => (e.getAttribute('aria-label') || '').includes('Scheme'))"))
        if not pos:
            raise RuntimeError("scheme switcher (role=combobox, aria-label ~ 'Scheme') not found")
        await cdp.click(*pos)
        for _poll in range(10):
            await asyncio.sleep(0.25)
            opt = await cdp.js(_center_js(
                "[...document.querySelectorAll('[role=\"option\"]')]"
                f".find(e => (e.textContent || '').trim() === {title_js})"))
            if opt:
                await cdp.click(*opt)
                return
        await cdp.send("Input.dispatchKeyEvent", type="keyDown", key="Escape")
    raise RuntimeError(f"scheme option {title!r} not found")


async def screenshot(cdp, path):
    res = await cdp.send("Page.captureScreenshot", format="png")
    path.write_bytes(base64.b64decode(res["data"]))


async def run_audit(url, out_dir, browser):
    inventory = scheme_inventory()
    proc, profile, port = launch_browser(browser)
    cdp = await CDP.connect(page_ws_url(port))
    results = {row["name"]: dict(row) for row in inventory}
    meta = {}
    try:
        await cdp.send("Page.enable")
        await cdp.send("Runtime.enable")

        # Desktop pass: one session, schemes in registry order.
        await cdp.send("Emulation.setDeviceMetricsOverride", **DESKTOP)
        await cdp.send("Page.navigate", url=url)
        first = inventory[0]["title"]
        results[inventory[0]["name"]]["ready_ms_desktop"] = await wait_ready(
            cdp, expect_title=first, min_wait=1.0)
        meta.update(await cdp.js(JS_TRANSFER))
        for idx, row in enumerate(inventory):
            if idx:
                await switch_scheme(cdp, row["title"])
                results[row["name"]]["ready_ms_desktop"] = await wait_ready(
                    cdp, expect_title=row["title"], min_wait=1.0)
            log(f"desktop: {row['title']}")
            results[row["name"]]["desktop"] = await cdp.js(JS_MEASURE)
            await screenshot(cdp, out_dir / f"{row['name']}_desktop.png")
            if row["n_advanced_controls"] and await cdp.js(JS_ADVANCED_TOGGLE):
                # The toggle may be client-side (expander) or a rerun (button):
                # wait for a settled page either way.
                await wait_ready(cdp, expect_title=row["title"], min_wait=1.2)
                m = await cdp.js(JS_MEASURE)
                results[row["name"]]["desktop_rail_advanced_open"] = m["rail_scroll_height"]
                await cdp.js(JS_ADVANCED_TOGGLE)
                await wait_ready(cdp, expect_title=row["title"], min_wait=0.8)

        # Mobile pass: fresh session at phone width, sidebar closed for measurement.
        await cdp.send("Emulation.setDeviceMetricsOverride", **MOBILE)
        await cdp.send("Page.navigate", url=url)
        await wait_ready(cdp, expect_title=first, min_wait=1.0)
        for idx, row in enumerate(inventory):
            # The switcher lives in the top bar; an open sidebar would cover it.
            await cdp.js(JS_SIDEBAR + "(false)")
            await asyncio.sleep(0.8)
            if idx:
                await switch_scheme(cdp, row["title"])
                await wait_ready(cdp, expect_title=row["title"], min_wait=1.0)
            log(f"mobile: {row['title']}")
            results[row["name"]]["mobile"] = await cdp.js(JS_MEASURE)
            await screenshot(cdp, out_dir / f"{row['name']}_mobile.png")
    except Exception:
        # Keep the page as it was at the failure for diagnosis.
        try:
            await screenshot(cdp, out_dir / "_failure.png")
        except Exception:
            pass
        raise
    finally:
        await cdp.close()
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except Exception:
            proc.kill()
        shutil.rmtree(profile, ignore_errors=True)
    return meta, results


def git_head():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=str(REPO),
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------
ROWS = [
    ("plot top (desktop, px)", lambda s: (s.get("desktop") or {}).get("plot", {}) and s["desktop"]["plot"]["top"]),
    ("plot fully on 1st screen", lambda s: (s.get("desktop") or {}).get("plot", {}) and s["desktop"]["plot"]["fully_visible"]),
    ("rail height (px)", lambda s: (s.get("desktop") or {}).get("rail_scroll_height")),
    ("rail content height", lambda s: (s.get("desktop") or {}).get("rail_content_height")),
    ("rail height, advanced open", lambda s: s.get("desktop_rail_advanced_open")),
    ("rail px / main control", lambda s: _per_control(s)),
    ("help glyphs", lambda s: (s.get("desktop") or {}).get("help_glyphs")),
    ("help ? icons", lambda s: (s.get("desktop") or {}).get("help_icons")),
    ("iframes", lambda s: (s.get("desktop") or {}).get("iframes")),
    ("main expanders", lambda s: (s.get("desktop") or {}).get("main_expanders")),
    ("text leaves, 1st screen", lambda s: (s.get("desktop") or {}).get("text_leaves_first_viewport")),
    ("plot top (phone, px)", lambda s: (s.get("mobile") or {}).get("plot", {}) and s["mobile"]["plot"]["top"]),
    ("ready after switch (ms)", lambda s: s.get("ready_ms_desktop")),
]


def _per_control(s):
    rail = (s.get("desktop") or {}).get("rail_scroll_height")
    n = s.get("n_main_controls")
    return round(rail / n) if rail and n else None


def summary_markdown(doc):
    names = list(doc["schemes"])
    head = "| metric | " + " | ".join(names) + " |"
    sep = "|---|" + "---|" * len(names)
    lines = [f"# UI audit `{doc['tag']}`", "",
             f"- commit `{doc['git_head']}` · {doc['created']} · {doc['url_mode']}",
             f"- initial transfer {doc['meta'].get('transfer_kb')} KB, {doc['meta'].get('requests')} requests",
             "- viewports: desktop 1440×900, phone 390×844", "", head, sep]
    for label, fn in ROWS:
        vals = []
        for n in names:
            try:
                v = fn(doc["schemes"][n])
            except Exception:
                v = None
            vals.append("—" if v is None or v == {} else str(v))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    return "\n".join(lines) + "\n"


def compare(tag_a, tag_b):
    a = json.loads((AUDIT_ROOT / tag_a / "metrics.json").read_text(encoding="utf-8"))
    b = json.loads((AUDIT_ROOT / tag_b / "metrics.json").read_text(encoding="utf-8"))
    names = [n for n in a["schemes"] if n in b["schemes"]]
    print(f"| metric | " + " | ".join(names) + " |")
    print("|---|" + "---|" * len(names))
    for label, fn in ROWS:
        cells = []
        for n in names:
            try:
                va, vb = fn(a["schemes"][n]), fn(b["schemes"][n])
            except Exception:
                va = vb = None
            cells.append(f"{va if va not in (None, {}) else '—'} → {vb if vb not in (None, {}) else '—'}")
        print(f"| {label} | " + " | ".join(cells) + " |")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tag", help="output folder under docs/ui_redesign/audit/")
    ap.add_argument("--url", help="audit an already running app instead of spawning one")
    ap.add_argument("--browser", help="Edge/Chrome executable")
    ap.add_argument("--compare", nargs=2, metavar=("TAG_A", "TAG_B"))
    args = ap.parse_args()
    if args.compare:
        compare(*args.compare)
        return
    if not args.tag:
        ap.error("--tag is required unless --compare is used")

    out_dir = AUDIT_ROOT / args.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    server = None
    url = args.url
    if url is None:
        server, url = spawn_streamlit()
    try:
        meta, results = asyncio.run(run_audit(url, out_dir, find_browser(args.browser)))
    finally:
        if server is not None:
            server.terminate()
            try:
                server.wait(timeout=10)
            except Exception:
                server.kill()
    doc = dict(tag=args.tag, git_head=git_head(),
               created=time.strftime("%Y-%m-%d %H:%M:%S %z"),
               url_mode="spawned" if server is not None else "external",
               meta=meta, schemes=results)
    (out_dir / "metrics.json").write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                                          encoding="utf-8")
    (out_dir / "summary.md").write_text(summary_markdown(doc), encoding="utf-8")
    print(summary_markdown(doc))


if __name__ == "__main__":
    main()
