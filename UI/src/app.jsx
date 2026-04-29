// Novum Pipeline — main app (single-scroll layout)

const { useState, useEffect, useRef, useMemo } = React;

function ActionsBar({ cfg, onSave, onRun, onToggleYAML, theme, onTheme }) {
  const queriesOk = cfg.queries.length >= 1 && cfg.queries.every((q) => q.key.trim() && q.path.trim());
  const targetsOk = cfg.targets.length >= 1 && cfg.targets.every((t) => t.key.trim() && t.path && t.path.trim());
  const canRun = queriesOk && targetsOk && cfg.results_dir.trim();
  return (
    <div className="actions-bar">
      <div className="status">
        <span className="status-dot" style={{ background: canRun ? "var(--accent)" : "var(--warn)" }} />
        {canRun ?
          <span>Ready · {cfg.queries.length} × query · {cfg.targets.length} × target · {cfg.rnahybrid.threads} threads</span> :
          <span>Resolve required fields to run</span>
        }
      </div>
      <button
        className="btn ghost sm"
        onClick={() => onTheme(theme === "dark" ? "light" : "dark")}
        title="Toggle theme"
        aria-label="Toggle theme"
        style={{ marginLeft: 10 }}
      >
        {theme === "dark" ? <SunIcon /> : <MoonIcon />}
      </button>
      <div style={{ flex: 1 }} />
      <div style={{ display: "flex", gap: 8 }}>
        <button className="btn" onClick={onToggleYAML}>
          <span style={{ width: 6, height: 6, borderRadius: "50%", background: "var(--accent)", boxShadow: "0 0 0 3px var(--accent-a18)", marginRight: 2 }} />
          <CodeIcon /> config.yaml
        </button>
        <button className="btn" onClick={onSave}><SaveIcon /> Save</button>
        <button className="btn primary" onClick={onRun} disabled={!canRun}>
          <PlayIcon /> Run pipeline
        </button>
      </div>
    </div>
  );
}

function YAMLDrawer({ open, onClose, cfg, onSave }) {
  const lines = useMemo(() => window.buildYAML(cfg), [cfg]);
  const html = useMemo(() => window.renderYAML(lines), [lines]);
  // `plain` is only needed when the user clicks Copy — compute lazily so the
  // drawer doesn't pay for it on every keystroke that mutates `cfg`.
  const copy = () => navigator.clipboard?.writeText(window.renderYAMLPlain(lines));

  return (
    <React.Fragment>
      <div className={`yaml-backdrop ${open ? "open" : ""}`} onClick={onClose} />
      <aside className={`yaml-drawer ${open ? "open" : ""}`} aria-hidden={!open}>
        <div className="head">
          <div>
            <div className="ttl">Generated config</div>
            <div className="path">Config/config.yaml</div>
          </div>
          <button className="btn ghost sm" onClick={onClose} aria-label="Close">✕</button>
        </div>
        <div className="body">
          <pre className="yaml-panel" dangerouslySetInnerHTML={{ __html: html }} />
        </div>
        <div className="foot">
          <button className="btn sm" onClick={copy}><CopyIcon /> Copy</button>
          <button className="btn primary sm" onClick={onSave}><SaveIcon /> Save to disk</button>
        </div>
      </aside>
    </React.Fragment>
  );
}

function Hero({ cfg }) {
  const setCount = window.RNAHYBRID_NULLABLE_KEYS.filter((k) => cfg.rnahybrid[k].set).length;
  return (
    <div style={{
      display: "flex", alignItems: "center", gap: 24,
      padding: "18px 22px",
      background: "linear-gradient(135deg, var(--accent-soft), transparent 70%), var(--surface)",
      border: "1px solid var(--line)",
      borderRadius: "var(--r-lg)",
      boxShadow: "var(--shadow-1)"
    }}>
      <div style={{ flex: 1 }}>
        <div style={{ fontSize: 11, letterSpacing: "0.08em", textTransform: "uppercase", color: "var(--accent-ink)", fontWeight: 600 }}>
          Configure run
        </div>
        <div style={{ fontSize: 20, fontWeight: 600, letterSpacing: "-0.01em", marginTop: 4 }}>
          RNA → bacterial target hybridization
        </div>
        <div style={{ fontSize: 13, color: "var(--fg-3)", marginTop: 4 }}>
          Build your pipeline config, then hit run.
        </div>
      </div>
      <div style={{ display: "flex", gap: 10 }}>
        <Stat label="Queries" value={cfg.queries.length} />
        <Stat label="Targets" value={cfg.targets.length} />
        <Stat label="Threads" value={cfg.rnahybrid.threads} />
        <Stat label="Set params" value={`${setCount}/${window.RNAHYBRID_NULLABLE_KEYS.length}`} />
      </div>
    </div>
  );
}

function Stat({ label, value }) {
  return (
    <div style={{
      padding: "10px 14px",
      background: "var(--surface)",
      border: "1px solid var(--line)",
      borderRadius: "var(--r-md)",
      minWidth: 76,
      textAlign: "center"
    }}>
      <div style={{ fontFamily: "var(--font-mono)", fontSize: 18, fontWeight: 600, color: "var(--fg)" }}>{value}</div>
      <div style={{ fontSize: 10.5, letterSpacing: "0.06em", textTransform: "uppercase", color: "var(--fg-4)", marginTop: 2 }}>{label}</div>
    </div>
  );
}

function App() {
  const [cfg, setCfg] = useState(() => JSON.parse(JSON.stringify(window.INITIAL_CONFIG)));
  const [theme, setTheme] = useState(() => localStorage.getItem("np:theme") || "light");
  const [yamlOpen, setYamlOpen] = useState(false);
  const [toast, setToast] = useState(null);
  const toastTimer = useRef(null);

  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    localStorage.setItem("np:theme", theme);
  }, [theme]);

  const showToast = (msg) => {
    clearTimeout(toastTimer.current);
    setToast(msg);
    toastTimer.current = setTimeout(() => setToast(null), 3200);
  };
  const serializeYAML = () => window.renderYAMLPlain(window.buildYAML(cfg));
  const save = async () => {
    const api = window.pywebview && window.pywebview.api;
    if (!api) return showToast("Preview mode · nothing written");
    const res = await api.save_config(serializeYAML());
    if (res && res.ok) showToast(`Saved to ${res.path}`);
    else showToast(`Save failed · ${res && res.error}`);
  };
  const run = async () => {
    const api = window.pywebview && window.pywebview.api;
    if (!api) return showToast("Preview mode · cannot launch");
    const res = await api.run_pipeline(serializeYAML());
    if (res && res.ok) showToast(`Pipeline launched · PID ${res.pid}`);
    else showToast(`Run failed · ${res && res.error}`);
  };

  return (
    <div className="desk">
      <div className="window" style={{ height: "calc(100vh / var(--ui-zoom))" }}>
        <div style={{ flex: 1, overflow: "hidden", display: "flex", flexDirection: "column", minHeight: 0, position: "relative" }}>
          <div className="no-scrollbar" style={{ flex: 1, overflowY: "auto" }}>
            <div style={{ padding: 24, display: "flex", flexDirection: "column", gap: 18, margin: "0 auto", width: "100%" }}>
              <Hero cfg={cfg} />
              <QueriesSection cfg={cfg} setCfg={setCfg} />
              <TargetsSection cfg={cfg} setCfg={setCfg} />
              <RNACalibrateSection cfg={cfg} setCfg={setCfg} />
              <RNAHybridSection cfg={cfg} setCfg={setCfg} />
              <BuildPlotsSection cfg={cfg} setCfg={setCfg} />
              <OutputSection cfg={cfg} setCfg={setCfg} />
            </div>
          </div>
          <ActionsBar
            cfg={cfg}
            onSave={save}
            onRun={run}
            onToggleYAML={() => setYamlOpen((o) => !o)}
            theme={theme}
            onTheme={setTheme}
          />
          <YAMLDrawer open={yamlOpen} onClose={() => setYamlOpen(false)} cfg={cfg} onSave={save} />
        </div>
      </div>

      {toast && (
        <div style={{
          position: "fixed", bottom: 96, right: 28,
          padding: "10px 14px", background: "var(--fg)", color: "var(--bg)",
          borderRadius: "var(--r-sm)", fontSize: 12.5, fontWeight: 500,
          boxShadow: "var(--shadow-pop)", zIndex: 50
        }}>
          {toast}
        </div>
      )}
    </div>
  );
}

const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
