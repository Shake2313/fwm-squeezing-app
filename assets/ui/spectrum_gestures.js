// Shared by spectrum views: preview projected coordinates on animation frames and
// commit once on pointer release. The caller owns the physics and what gets painted.
export function bindSpectrumGesture(surface, state, {begin, preview, commit, cancel}) {
  let frame = 0, pending = null;
  function clearFrame() {
    if (frame) cancelAnimationFrame(frame);
    frame = 0;
    pending = null;
  }
  function start(event) {
    if (event.button !== 0 || event.isPrimary === false || state.drag) return;
    const drag = begin(event);
    if (!drag) return;
    // The pointer moves chart coordinates, not the browser's text selection.
    // Do this only after hit-testing so nearby editable fields work normally.
    event.preventDefault?.();
    state.drag = {...drag, pointerId:event.pointerId};
    surface.setPointerCapture?.(event.pointerId);
    preview(state.drag, event);
  }
  function move(event) {
    if (!state.drag || event.pointerId !== state.drag.pointerId) return;
    pending = event;
    if (!frame) frame = requestAnimationFrame(() => {
      frame = 0;
      if (state.drag && pending) preview(state.drag, pending);
      pending = null;
    });
  }
  function finish(event, cancelled) {
    if (!state.drag || event.pointerId !== state.drag.pointerId) return;
    const drag = state.drag;
    clearFrame();
    if (!cancelled) preview(drag, event);
    state.drag = null;
    if (surface.hasPointerCapture?.(event.pointerId)) surface.releasePointerCapture(event.pointerId);
    if (cancelled) cancel(drag);
    else commit(drag, event);
  }
  const up = event => finish(event, false), abort = event => finish(event, true);
  const escape = event => {
    if (event.key === 'Escape' && state.drag) {
      event.preventDefault();
      finish({pointerId:state.drag.pointerId}, true);
    }
  };
  surface.onpointerdown = start;
  // The chart may receive new HTML while a gesture is active. Its owner keeps
  // state.drag so new listeners can finish a capture from the old wrapper.
  window.addEventListener('pointermove', move, true);
  window.addEventListener('pointerup', up, true);
  window.addEventListener('pointercancel', abort, true);
  window.addEventListener('keydown', escape, true);
  return () => {
    clearFrame();
    if (surface.onpointerdown === start) surface.onpointerdown = null;
    window.removeEventListener('pointermove', move, true);
    window.removeEventListener('pointerup', up, true);
    window.removeEventListener('pointercancel', abort, true);
    window.removeEventListener('keydown', escape, true);
  };
}

export function spectrumPointer(svg, event, width, height) {
  const rect = svg.getBoundingClientRect();
  return {x:(event.clientX - rect.left) * width / rect.width,
          y:(event.clientY - rect.top) * height / rect.height};
}

const clampSpectrum = (value, low, high) => Math.min(high, Math.max(low, value));

// Both edges and the read cursor move in physical frequency units. Translation
// is clamped as a unit, so reaching a scan boundary never narrows the window.
export function navigateSpectrumWindow(action, initial, frequency, domain, minWidth) {
  const [a,b] = initial.window, [lo,hi] = domain, width = b-a;
  const x = clampSpectrum(frequency,lo,hi), shift = x-initial.x;
  let low=a, high=b, cursor=initial.cursor;
  switch (action) {
    case 'pan': {
      const delta = clampSpectrum(shift, lo-a, hi-b);
      low=a+delta; high=b+delta; cursor+=delta;
      break;
    }
    case 'center':
      low = clampSpectrum(x - width/2, lo, hi-width);
      high = low + width; cursor=x;
      break;
    case 'left': low=clampSpectrum(x,lo,b-minWidth); break;
    case 'right': high=clampSpectrum(x,a+minWidth,hi); break;
    case 'cursor': cursor=x; break;
  }
  low=clampSpectrum(low,lo,hi); high=clampSpectrum(high,lo,hi);
  return {window:[low,high], cursor:clampSpectrum(cursor,low,high)};
}
