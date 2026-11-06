// Rendering status only: no numerical values or Streamlit state are changed.
export function updateMessage(running, localPhase, failed=false) {
  if (failed) return {phase:'error', text:'Update failed · previous results may still be shown.'};
  if (localPhase==='preview') return {phase:'preview', text:'Preview · result cards show the last committed values.'};
  if (running || localPhase==='syncing') return {phase:'syncing', text:'Updating results… previous values may still be shown.'};
  // Keep the reserved status slot, but a synchronized result needs no notice.
  return {phase:'ready', text:''};
}

export default function({parentElement, data}) {
  const root=parentElement.querySelector('.update-status');
  if(!root) return;
  const app=document.querySelector('[data-testid="stApp"]');
  const state=root._updateState || (root._updateState={});
  state.cleanup?.();
  if(state.scheme!==data.scheme) {state.localPhase='ready';state.plotPending=false;state.plotOwner=null;}
  state.scheme=data.scheme;
  const paint=()=> {
    const running=app?.getAttribute('data-test-script-state')==='running';
    const failed=!!document.querySelector('[data-testid="stException"]');
    const message=updateMessage(running,state.localPhase,failed);
    if(root.textContent!==message.text) root.textContent=message.text;
    root.dataset.phase=message.phase;
  };
  const onPlot=e=> {
    if(e.detail?.scheme!==data.scheme) return;
    state.plotOwner=e.detail.owner;
    state.plotPending=e.detail.phase!=='ready';
    state.localPhase=e.detail.phase;
    paint();
  };
  const onInput=()=> {state.localPhase='syncing';paint();};
  let wasRunning=app?.getAttribute('data-test-script-state')==='running';
  // Keep the native Streamlit image and its full-screen control. The plot and
  // its description can arrive on separate deltas during a fragment rerun.
  const labelNativePlot=()=> {
    const card=app?.querySelector?.('.st-key-gabes_plotcard');
    if(!card) return;
    const description=card.querySelector('[data-spectrum-alt]')?.dataset.spectrumAlt;
    const image=card.querySelector('.st-key-gabes_plotstage [data-testid="stImage"] img');
    if(image && description && image.alt!==description) image.alt=description;
  };
  const observer=new MutationObserver((records)=> {
    if(records?.length && records.every(record=>record.type==='childList')) {
      if(records.some(record=>record.target.closest?.('.st-key-gabes_plotcard'))) labelNativePlot();
      return;
    }
    labelNativePlot();
    const running=app?.getAttribute('data-test-script-state')==='running';
    // A graph owns its preview/queue until its own ready event. Ordinary
    // controls are acknowledged when Streamlit completes the script.
    if(wasRunning && !running) {
      // A mode switch may remove the live plot without producing a ready
      // event. A detached owner must not keep this page in Updating forever.
      if(state.plotOwner && !state.plotOwner.isConnected) state.plotPending=false;
      if(!state.plotPending) state.localPhase='ready';
    }
    wasRunning=running;
    paint();
  });
  window.addEventListener('gabes-plot-status',onPlot);
  window.addEventListener('gabes-input-commit',onInput);
  if(app) observer.observe(app,{attributes:true,attributeFilter:['data-test-script-state','data-spectrum-alt'],childList:true,subtree:true});
  labelNativePlot();
  paint();
  const cleanup=()=> {
    observer.disconnect();
    window.removeEventListener('gabes-plot-status',onPlot);
    window.removeEventListener('gabes-input-commit',onInput);
  };
  state.cleanup=cleanup;
  return cleanup;
}
