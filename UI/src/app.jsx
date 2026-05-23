// Novum Pipeline — main app (single-scroll layout)

const { useState, useEffect, useRef } = React;

// Pipeline lifecycle:
//   idle      → no run in flight
//   running   → snakemake is alive
//   succeeded → terminal · exit 0
//   failed    → terminal · non-zero exit
//   cancelled → terminal · user pressed Cancel
//
// Terminal states stay until the user dismisses the overlay; we then call
// `acknowledge_pipeline` so the backend resets to idle. Without that
// handshake a fast finish could be missed between two polls.

function App() {
  const [cfg, setCfg] = useState(() => JSON.parse(JSON.stringify(window.INITIAL_CONFIG)));
  const [theme, setTheme] = useState(() => localStorage.getItem("np:theme") || "light");
  const [yamlOpen, setYamlOpen] = useState(false);
  const [toast, setToast] = useState(null);
  const toastTimer = useRef(null);

  const [savingState, setSavingState] = useState("idle");
  const [pipelineStatus, setPipelineStatus] = useState(null);
  const [overlayOpen, setOverlayOpen] = useState(false);

  const pipelineState = pipelineStatus?.state ?? "idle";

  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    localStorage.setItem("np:theme", theme);
  }, [theme]);

  useEffect(() => () => clearTimeout(toastTimer.current), []);

  const showToast = (msg, tone = "neutral") => {
    clearTimeout(toastTimer.current);
    setToast({ msg, tone });
    toastTimer.current = setTimeout(() => setToast(null), 3200);
  };

  const serializeYAML = () => window.renderYAMLPlain(window.buildYAML(cfg));

  const save = async () => {
    if (savingState === "saving") return;
    const api = window.pywebview?.api;
    setSavingState("saving");
    try {
      if (!api) {
        showToast("Preview mode · nothing written");
        return;
      }
      const res = await api.save_config(serializeYAML());
      if (res?.ok) showToast(`Saved to ${res.path}`, "success");
      else showToast(`Save failed · ${res?.error || "unknown error"}`, "danger");
    } finally {
      setSavingState("idle");
    }
  };

  const run = async () => {
    if (pipelineState === "running") return;
    const api = window.pywebview?.api;
    if (!api?.run_pipeline) {
      showToast("Preview mode · cannot launch");
      return;
    }
    // Clear any stale terminal state so the polling effect can't snap the
    // overlay back to a prior 'succeeded'/'failed' as the new run starts.
    // If the backend refuses (a previous run is still draining), surface
    // that instead of plowing into run_pipeline which would also reject.
    if (api.acknowledge_pipeline) {
      const ack = await api.acknowledge_pipeline();
      if (ack && ack.ok === false) {
        showToast("Previous run still finishing — try again", "danger");
        return;
      }
    }
    setOverlayOpen(true);
    setPipelineStatus({ state: "running", elapsed: 0 });
    const res = await api.run_pipeline(serializeYAML());
    if (!res?.ok) {
      setPipelineStatus({ state: "failed", returncode: -1, elapsed: 0 });
      showToast(`Run failed · ${res?.error || "unknown error"}`, "danger");
    }
  };

  const cancel = async () => {
    const api = window.pywebview?.api;
    if (api?.cancel_pipeline) await api.cancel_pipeline();
  };

  const dismissOverlay = async () => {
    const api = window.pywebview?.api;
    if (api?.acknowledge_pipeline) await api.acknowledge_pipeline();
    setOverlayOpen(false);
    setPipelineStatus(null);
  };

  useEffect(() => {
    if (pipelineState !== "running") return;
    const api = window.pywebview?.api;
    if (!api?.pipeline_status) return;
    let cancelled = false;

    const poll = async () => {
      if (cancelled) return;
      try {
        const status = await api.pipeline_status();
        // Race: between optimistic running and run_pipeline returning, the
        // backend reports idle. Ignore so the overlay doesn't flicker back.
        if (cancelled || status.state === "idle") return;
        setPipelineStatus(status);
        if (status.state === "succeeded") {
          showToast(`Pipeline complete · ${window.formatElapsed(status.elapsed)}`, "success");
        } else if (status.state === "failed") {
          showToast(`Pipeline failed · exit ${status.returncode}`, "danger");
        } else if (status.state === "cancelled") {
          showToast("Pipeline cancelled");
        }
      } catch (_) { /* transient bridge errors retry on the next tick */ }
    };

    poll();
    const id = setInterval(poll, 1000);
    return () => { cancelled = true; clearInterval(id); };
  }, [pipelineState]);

  return (
    <div className="desk">
      <div className="window" style={{ height: "calc(100vh / var(--ui-zoom))" }}>
        <div style={{ flex: 1, overflow: "hidden", display: "flex", flexDirection: "column", minHeight: 0, position: "relative" }}>
          <div className="no-scrollbar" style={{ flex: 1, overflowY: "auto" }}>
            <div style={{ padding: 24, display: "flex", flexDirection: "column", gap: 18, margin: "0 auto", width: "100%" }}>
              <ErrorBoundary>
                <Hero cfg={cfg} />
                <QueriesSection cfg={cfg} setCfg={setCfg} />
                <TargetsSection cfg={cfg} setCfg={setCfg} />
                <SharedParamsSection cfg={cfg} setCfg={setCfg} />
                <RNACalibrateSection cfg={cfg} setCfg={setCfg} />
                <RNAHybridSection cfg={cfg} setCfg={setCfg} />
                <IntaRNASection cfg={cfg} setCfg={setCfg} />
                <PlotsSection cfg={cfg} setCfg={setCfg} />
                <OutputSection cfg={cfg} setCfg={setCfg} />
              </ErrorBoundary>
            </div>
          </div>
          <ActionsBar
            cfg={cfg}
            onSave={save}
            onRun={run}
            onToggleYAML={() => setYamlOpen((o) => !o)}
            theme={theme}
            onTheme={setTheme}
            pipelineState={pipelineState}
            elapsed={pipelineStatus?.elapsed}
            savingState={savingState}
          />
          <YAMLDrawer open={yamlOpen} onClose={() => setYamlOpen(false)} cfg={cfg} onSave={save} />
          <RunOverlay
            open={overlayOpen}
            status={pipelineStatus}
            threads={cfg.threads}
            queries={cfg.queries.length}
            targets={cfg.targets.length}
            onCancel={cancel}
            onClose={dismissOverlay}
          />
        </div>
      </div>

      {toast && <Toast tone={toast.tone}>{toast.msg}</Toast>}
    </div>
  );
}

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
