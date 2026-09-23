"""User's Guide launcher — an inline `st.components.v2` button (no iframe).

Streamlit's static handler serves the guide HTML as text/plain, so a plain
link would show source. The button fetches the file and opens it as a
text/html Blob in a new tab (the 851 KB guide is fetched only on click).
"""
import streamlit as st

_HTML = """<button type="button" class="gbtn" aria-label="User's Guide" title="사용자 안내서를 새 창에서 엽니다">
<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linejoin="round" aria-hidden="true"><path d="M8 4.5C6.8 3.6 5 3.2 2.5 3.3v8.4c2.5-.1 4.3.3 5.5 1.2 1.2-.9 3-1.3 5.5-1.2V3.3C11 3.2 9.2 3.6 8 4.5z"/><path d="M8 4.5v8.4"/></svg>
<span class="lbl">Guide</span></button>"""

_CSS = """
.gbtn{display:inline-flex;align-items:center;gap:.4rem;height:2rem;padding:0 .65rem;
  border:1px solid transparent;border-radius:6px;background:transparent;cursor:pointer;
  color:var(--g-ink-2,#1F2D3D);font:500 .84rem/1 var(--g-font-sans,sans-serif);}
.gbtn:hover{background:var(--g-fill-soft,#EDF1F5);}
.gbtn:focus-visible{outline:2px solid var(--g-accent,#0284C7);outline-offset:2px;}
.gbtn svg{width:16px;height:16px;flex-shrink:0;}
.gbtn[aria-busy="true"]{opacity:.6;cursor:progress;}
@media (max-width:640px){.gbtn .lbl{display:none;}}
"""

_JS = """
export default function(component) {
  const { parentElement, data } = component;
  const btn = parentElement.querySelector("button");
  if (!btn) return;
  const url = (data && data.url) || "app/static/GABES_User_Guide.html";
  btn.onclick = async () => {
    btn.setAttribute("aria-busy", "true");
    try {
      const text = await (await fetch(url)).text();
      const blobUrl = URL.createObjectURL(new Blob([text], { type: "text/html" }));
      if (!window.open(blobUrl, "_blank")) window.open(url, "_blank");
    } catch (_e) {
      window.open(url, "_blank");
    } finally {
      btn.removeAttribute("aria-busy");
    }
  };
}
"""

_GUIDE_BUTTON = st.components.v2.component("gabes_guide_button", html=_HTML, css=_CSS, js=_JS)


def guide_button(url, key="gabes_guide"):
    """Mount the launcher in the current container."""
    return _GUIDE_BUTTON(key=key, data={"url": url}, width="content")
