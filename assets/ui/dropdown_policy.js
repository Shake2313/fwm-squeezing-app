// Native Streamlit selectboxes are choice controls. Keep their keyboard
// navigation while closing the text-search editor for these short option lists.
const dropdownInput = '[data-testid="stSelectbox"] input[role="combobox"]';

export default function mountDropdownPolicy({parentElement}) {
  const marker = parentElement.querySelector('.dropdown-policy');
  if (!marker) return;
  marker._cleanup?.();

  function prepare() {
    for (const input of document.querySelectorAll(dropdownInput)) {
      input.readOnly = true;
      input.setAttribute('inputmode', 'none');
      input.setAttribute('aria-autocomplete', 'none');
    }
  }
  function target(event) {
    return event.target?.matches?.(dropdownInput);
  }
  function stopText(event) {
    if (!target(event)) return;
    event.preventDefault();
    event.stopPropagation();
  }
  function keydown(event) {
    if (!target(event)) return;
    const textKey = event.key.length === 1 && event.key !== ' ';
    const deletion = event.key === 'Backspace' || event.key === 'Delete';
    // Preserve shortcuts such as Tab, Escape and copy, but never paste/search.
    const paste = (event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'v';
    if ((!event.ctrlKey && !event.metaKey && !event.altKey && textKey) || deletion || paste) {
      stopText(event);
    }
  }
  const observer = new MutationObserver(prepare);
  observer.observe(document.body, {childList:true, subtree:true});
  document.addEventListener('keydown', keydown, true);
  for (const type of ['beforeinput', 'paste', 'drop']) {
    document.addEventListener(type, stopText, true);
  }
  prepare();
  const cleanup = () => {
    observer.disconnect();
    document.removeEventListener('keydown', keydown, true);
    for (const type of ['beforeinput', 'paste', 'drop']) {
      document.removeEventListener(type, stopText, true);
    }
  };
  marker._cleanup = cleanup;
  return cleanup;
}
