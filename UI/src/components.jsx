// Novum Pipeline — shared form components (React, Babel)

const { useState, useRef, useEffect, useMemo } = React;

// --- primitives ---------------------------------------------------

function Toggle({ on, onChange, ariaLabel, disabled = false }) {
  return (
    <button
      type="button"
      className={`toggle ${on ? "on" : ""}`}
      onClick={() => !disabled && onChange(!on)}
      disabled={disabled}
      aria-label={ariaLabel}
      aria-pressed={on}
    />
  );
}

// Grey labeled subcard. `right` lets a Toggle/badge sit on the label row,
// mirroring NullableField's top-row treatment.
function Subcard({ label, hint, right, children }) {
  return (
    <div className="subcard">
      <div className="subcard-label">
        <span>{label}{hint && <em> · {hint}</em>}</span>
        {right}
      </div>
      {children}
    </div>
  );
}

// Subcard whose body is a flex-wrap row of PillToggles or similar pills.
// Used for the three IntaRNA "Constraints" blocks and any future pill row.
function PillGroup({ label, hint, children }) {
  return (
    <Subcard label={label} hint={hint}>
      <div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>
        {children}
      </div>
    </Subcard>
  );
}

// Pill-style boolean toggle — same visual language as the Plots "Plot types"
// selector. Used for IntaRNA constraint groups (seed / accessibility / output).
function PillToggle({ on, onChange, label }) {
  return (
    <button
      type="button"
      onClick={() => onChange(!on)}
      className="btn sm"
      aria-pressed={on}
      style={{
        background: on ? "var(--accent-soft)" : "var(--surface)",
        color: on ? "var(--accent-ink)" : "var(--fg-2)",
        borderColor: on ? "var(--accent)" : "var(--line)",
      }}
    >
      {on ? "✓ " : ""}{label}
    </button>
  );
}

// Single-select segmented control. A radiogroup (not a tablist: there are no
// tabpanels) with roving tabindex and arrow-key navigation, so it's one tab
// stop and screen readers announce the selected option correctly.
function SegmentedControl({ options, value, onChange, ariaLabel }) {
  const ref = useRef(null);
  const idx = Math.max(0, options.findIndex((o) => o.value === value));
  const onKeyDown = (e) => {
    const dir = (e.key === "ArrowRight" || e.key === "ArrowDown") ? 1
              : (e.key === "ArrowLeft" || e.key === "ArrowUp") ? -1 : 0;
    if (!dir) return;
    e.preventDefault();
    const next = (idx + dir + options.length) % options.length;
    onChange(options[next].value);
    ref.current?.querySelectorAll('[role="radio"]')[next]?.focus({ preventScroll: true });
  };
  return (
    <div className="variation-switch" role="radiogroup" aria-label={ariaLabel} ref={ref} onKeyDown={onKeyDown}>
      {options.map((o) => {
        const active = o.value === value;
        return (
          <button
            key={o.value}
            type="button"
            role="radio"
            aria-checked={active}
            tabIndex={active ? 0 : -1}
            className={active ? "active" : ""}
            onClick={() => onChange(o.value)}
            style={{ padding: "6px 14px" }}
          >
            {o.label}
          </button>
        );
      })}
    </div>
  );
}

// Focus management for the run overlay and YAML drawer. On open: remember the
// previously focused element, move focus inside, trap Tab, and (when onEscape
// is supplied) close on Escape. On close: restore focus to where it was. Esc
// behaviour is read through a ref so callers can gate it (e.g. no dismiss while
// the pipeline is still running) without re-running the trap.
function useDialog(ref, open, { onEscape } = {}) {
  const prevFocus = useRef(null);
  const escRef = useRef(onEscape);
  escRef.current = onEscape;

  useEffect(() => {
    if (!open) return;
    const node = ref.current;
    if (!node) return;
    prevFocus.current = document.activeElement;
    const SEL = 'a[href],button:not([disabled]),input:not([disabled]),select:not([disabled]),textarea:not([disabled]),[tabindex]:not([tabindex="-1"])';
    const focusables = () => Array.from(node.querySelectorAll(SEL));
    // preventScroll: the drawer/overlay is absolutely positioned and animates
    // in via transform; a plain focus() would scroll the (overflow:hidden but
    // still scrollable) content area to "reach" it, jolting the background.
    (focusables()[0] || node).focus({ preventScroll: true });

    const onKeyDown = (e) => {
      if (e.key === "Escape") {
        if (escRef.current) { e.stopPropagation(); escRef.current(); }
        return;
      }
      if (e.key !== "Tab") return;
      const items = focusables();
      if (items.length === 0) { e.preventDefault(); return; }
      const first = items[0], last = items[items.length - 1];
      if (e.shiftKey && document.activeElement === first) {
        e.preventDefault(); last.focus({ preventScroll: true });
      } else if (!e.shiftKey && document.activeElement === last) {
        e.preventDefault(); first.focus({ preventScroll: true });
      }
    };
    node.addEventListener("keydown", onKeyDown);
    return () => {
      node.removeEventListener("keydown", onKeyDown);
      const prev = prevFocus.current;
      prevFocus.current = null;
      if (prev && prev.focus && document.contains(prev)) prev.focus({ preventScroll: true });
    };
  }, [open]);
}

// Catches render errors anywhere in its subtree. Without this, a thrown
// exception unmounts the whole App and the user sees a blank window.
class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { error: null, info: null };
  }
  static getDerivedStateFromError(error) { return { error }; }
  componentDidCatch(error, info) {
    this.setState({ error, info });
    console.error("UI error:", error, info);
  }
  render() {
    if (!this.state.error) return this.props.children;
    return (
      <div style={{
        margin: 24, padding: 18,
        background: "var(--danger-a12)",
        border: "1px solid var(--danger)",
        borderRadius: "var(--r-sm)",
        fontFamily: "var(--font-mono)",
        fontSize: 12.5,
        color: "var(--fg)",
        whiteSpace: "pre-wrap",
      }}>
        <div style={{ fontWeight: 700, marginBottom: 8 }}>UI render error</div>
        <div>{String(this.state.error)}</div>
        {this.state.info?.componentStack && (
          <div style={{ color: "var(--fg-3)", marginTop: 8 }}>
            {this.state.info.componentStack}
          </div>
        )}
        <button className="btn sm" style={{ marginTop: 12 }}
          onClick={() => this.setState({ error: null, info: null })}>
          Reset
        </button>
      </div>
    );
  }
}

// Browser-preview fallback when no native file picker is available — used by
// PathInput (single path) and KeyedFileRow's free-form mode.
function promptBrowse(current, placeholder, label = "Select…") {
  return prompt(`${label} (simulated)`, current || placeholder || "");
}

function PathInput({ value, onChange, placeholder, mono = true }) {
  return (
    <div className="field-row">
      <input
        className={`input ${mono ? "mono" : ""}`}
        value={value}
        placeholder={placeholder}
        onChange={(e) => onChange(e.target.value)}
        spellCheck={false}
      />
      <button
        className="btn"
        type="button"
        onClick={() => {
          const pick = promptBrowse(value, placeholder, "Select folder");
          if (pick != null) onChange(pick);
        }}
        title="Browse…"
      >
        <FolderIcon /> Browse…
      </button>
    </div>
  );
}

// Plain helper (no React state) — generates the add/update/remove closures
// every list-edit section needs. `minLength` keeps a section from emptying
// itself below a sentinel (queries/targets require at least one row).
function useListEditor(items, setItems, makeItem, minLength = 0) {
  return {
    add: () => setItems([...items, makeItem()]),
    update: (idx, patch) => {
      const copy = items.slice();
      copy[idx] = { ...copy[idx], ...patch };
      setItems(copy);
    },
    remove: (idx) => {
      if (items.length <= minLength) return;
      const copy = items.slice();
      copy.splice(idx, 1);
      setItems(copy);
    },
  };
}

// --- launcher bridge ---------------------------------------------

// Calls a Python API method exposed by UI/launcher.py (via PyWebView's js_api).
// Returns { items, hasApi }:
//   items === null  → still loading (waiting for pywebviewready)
//   items === []    → API returned an empty list, OR preview mode (no API)
//   items === [...] → API returned files
//   hasApi          → true when window.pywebview.api[method] exists
function useApiList(method) {
  const [items, setItems] = useState(null);
  const [hasApi, setHasApi] = useState(false);

  useEffect(() => {
    let cancelled = false;
    const fetchOnce = () => {
      const api = window.pywebview && window.pywebview.api;
      if (!api || typeof api[method] !== "function") return false;
      setHasApi(true);
      Promise.resolve(api[method]()).then((data) => {
        if (!cancelled) setItems(Array.isArray(data) ? data : []);
      }).catch(() => {
        if (!cancelled) setItems([]);
      });
      return true;
    };

    if (fetchOnce()) {
      return () => { cancelled = true; };
    }

    const onReady = () => fetchOnce();
    window.addEventListener("pywebviewready", onReady);

    // Browser-preview safety net: if the bridge never appears, fall back.
    const timeoutId = setTimeout(() => {
      if (!cancelled && !window.pywebview) {
        setItems([]);
        setHasApi(false);
      }
    }, 600);

    return () => {
      cancelled = true;
      clearTimeout(timeoutId);
      window.removeEventListener("pywebviewready", onReady);
    };
  }, [method]);

  return { items, hasApi };
}

// --- keyed file/path row -----------------------------------------

// One row used by both Queries and Targets. When `items` is provided
// (an array, possibly empty, or null while loading) it renders the
// launcher-driven dropdown picker; otherwise it falls back to a free-form
// path input + simulated browse for browser previews.
function KeyedFileRow({
  keyVal, pathVal, onKey, onPath, onRemove, canRemove,
  placeholderKey, placeholderPath, items, emptyHint,
}) {
  const usePicker = items !== undefined;
  const cols = usePicker ? "180px 1fr auto" : "180px 1fr auto auto";
  const loading = items === null;
  const empty = usePicker && !loading && items.length === 0;

  return (
    <div style={{ display: "grid", gridTemplateColumns: cols, gap: 8, alignItems: "center" }}>
      <input
        className="input mono"
        value={keyVal}
        placeholder={placeholderKey}
        onChange={(e) => onKey(e.target.value)}
        aria-label="Key"
      />
      {usePicker ? (
        <select
          className="select"
          value={pathVal || ""}
          onChange={(e) => onPath(e.target.value)}
          disabled={loading || empty}
          aria-label="File"
        >
          <option value="" disabled>
            {loading ? "Loading…" : empty ? (emptyHint || "No files") : "— select a file —"}
          </option>
          {!loading && !empty && items.map((f) => (
            <option key={f.path} value={f.path}>{f.name}</option>
          ))}
        </select>
      ) : (
        <input
          className="input mono"
          value={pathVal}
          placeholder={placeholderPath}
          onChange={(e) => onPath(e.target.value)}
          spellCheck={false}
          aria-label="Path"
        />
      )}
      {!usePicker && (
        <button
          className="btn sm"
          type="button"
          onClick={() => {
            const pick = promptBrowse(pathVal, placeholderPath, "Select file");
            if (pick != null) onPath(pick);
          }}
          title="Browse…"
          aria-label="Browse…"
        >
          <FolderIcon />
        </button>
      )}
      <button
        className="btn sm danger-ghost"
        onClick={onRemove}
        disabled={!canRemove}
        title="Remove"
        aria-label="Remove"
      >
        <TrashIcon />
      </button>
    </div>
  );
}

// --- icons --------------------------------------------------------

function FolderIcon() {
  return (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden>
      <path d="M1.5 4.5a1.5 1.5 0 0 1 1.5-1.5h3l1.5 1.5h5a1.5 1.5 0 0 1 1.5 1.5V12a1.5 1.5 0 0 1-1.5 1.5H3A1.5 1.5 0 0 1 1.5 12V4.5Z"
        stroke="currentColor" strokeWidth="1.3" strokeLinejoin="round"/>
    </svg>
  );
}
function PlusIcon() { return <svg width="12" height="12" viewBox="0 0 12 12" aria-hidden><path d="M6 2v8M2 6h8" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round"/></svg>; }
function TrashIcon() { return <svg width="13" height="13" viewBox="0 0 16 16" fill="none" aria-hidden><path d="M3 5h10M6.5 5V3.5A.5.5 0 0 1 7 3h2a.5.5 0 0 1 .5.5V5M4.5 5l.5 8a1 1 0 0 0 1 1h4a1 1 0 0 0 1-1l.5-8" stroke="currentColor" strokeWidth="1.3" strokeLinecap="round" strokeLinejoin="round"/></svg>; }
function CodeIcon() { return <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden><path d="M6 4.5L2.5 8 6 11.5M10 4.5L13.5 8 10 11.5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/></svg>; }
function PlayIcon() { return <svg width="12" height="12" viewBox="0 0 12 12" aria-hidden><path d="M3 2.5v7l6-3.5-6-3.5Z" fill="currentColor"/></svg>; }
function SaveIcon() { return <svg width="13" height="13" viewBox="0 0 16 16" fill="none" aria-hidden><path d="M3 3h8.5L13 4.5V13a.5.5 0 0 1-.5.5h-9A.5.5 0 0 1 3 13V3Z" stroke="currentColor" strokeWidth="1.3" strokeLinejoin="round"/><path d="M5.5 3v3h4V3M5.5 13v-4h4v4" stroke="currentColor" strokeWidth="1.3" strokeLinejoin="round"/></svg>; }
function SunIcon() { return <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden><circle cx="8" cy="8" r="3" stroke="currentColor" strokeWidth="1.3"/><path d="M8 1.5v1.5M8 13v1.5M1.5 8H3M13 8h1.5M3.3 3.3l1 1M11.7 11.7l1 1M3.3 12.7l1-1M11.7 4.3l1-1" stroke="currentColor" strokeWidth="1.3" strokeLinecap="round"/></svg>; }
function MoonIcon() { return <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden><path d="M13 9.5A5.5 5.5 0 0 1 6.5 3c0-.5.07-1 .2-1.5A6 6 0 1 0 14.5 9.3c-.5.13-1 .2-1.5.2Z" fill="currentColor"/></svg>; }
function CopyIcon() { return <svg width="13" height="13" viewBox="0 0 16 16" fill="none" aria-hidden><rect x="4" y="4" width="9" height="9" rx="1.5" stroke="currentColor" strokeWidth="1.3"/><path d="M3 10.5V3.5A.5.5 0 0 1 3.5 3h7" stroke="currentColor" strokeWidth="1.3" strokeLinecap="round"/></svg>; }
function CheckIcon({ size = 14 }) { return <svg width={size} height={size} viewBox="0 0 16 16" fill="none" aria-hidden><path d="M3 8.5l3 3 7-7" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round"/></svg>; }
function XIcon({ size = 14 }) { return <svg width={size} height={size} viewBox="0 0 16 16" fill="none" aria-hidden><path d="M4 4l8 8M12 4l-8 8" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round"/></svg>; }
function StopIcon() { return <svg width="11" height="11" viewBox="0 0 12 12" aria-hidden><rect x="2.5" y="2.5" width="7" height="7" rx="1" fill="currentColor"/></svg>; }

// --- nullable numeric input --------------------------------------

function NullableField({
  name, doc, value, onChange,
  kind = "number",        // "number" | "integer" | "pair"
  placeholder,
  pairLabels,             // ["start nt", "end nt"] etc.
  forcedNull = false,
  forcedNullHint,
  nullHint,
  min,
  max
}) {
  const on = !forcedNull && value.set;
  // `<input type=number>` doesn't clamp typed input — bounce negatives only,
  // don't bump partial digits up to `min` mid-stroke.
  const clamp = (n) =>
    typeof n === "number" && min != null && min >= 0 && n < 0 ? min : n;
  const toInt = (s) => {
    if (s === "" || s === "-") return s;
    const n = parseInt(s, 10);
    return Number.isFinite(n) ? clamp(n) : "";
  };
  return (
    <div className={`nullable ${on ? "on" : ""} ${forcedNull ? "locked" : ""}`}>
      <div className="top">
        <div className="name">
          {name}{doc && <em> · {doc}</em>}
        </div>
        <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
          <span className="auto-tag">
            {forcedNull ? "LOCKED · null" : (on ? "SET" : "AUTO (null)")}
          </span>
          <Toggle
            on={on}
            onChange={(v) => onChange({ ...value, set: v })}
            disabled={forcedNull}
            ariaLabel={`Enable ${name}`}
          />
        </div>
      </div>
      {forcedNull ? (
        <div className="ghost-input">{forcedNullHint || "null · provided elsewhere"}</div>
      ) : on ? (
        kind === "pair" ? (
          <div className="pair-input">
            <input
              className="input mono"
              type="number"
              step="1"
              inputMode="numeric"
              min={min}
              max={max}
              placeholder={pairLabels?.[0]}
              value={value.a ?? ""}
              onChange={(e) => onChange({ ...value, a: toInt(e.target.value) })}
            />
            <span className="pair-sep">,</span>
            <input
              className="input mono"
              type="number"
              step="1"
              inputMode="numeric"
              min={min}
              max={max}
              placeholder={pairLabels?.[1]}
              value={value.b ?? ""}
              onChange={(e) => onChange({ ...value, b: toInt(e.target.value) })}
            />
          </div>
        ) : (
          <input
            className="input mono"
            type="number"
            inputMode={kind === "integer" ? "numeric" : "decimal"}
            step={kind === "integer" ? "1" : "any"}
            min={min}
            max={max}
            value={value.value}
            placeholder={placeholder}
            onChange={(e) => {
              const raw = e.target.value;
              if (raw === "") return onChange({ ...value, value: "" });
              if (kind === "integer") return onChange({ ...value, value: toInt(raw) });
              const n = Number(raw);
              // Ignore partial input like "1e" or ".": NaN would otherwise
              // serialize as the literal `NaN` in YAML.
              if (!Number.isFinite(n)) return;
              onChange({ ...value, value: clamp(n) });
            }}
          />
        )
      ) : (
        <div className="ghost-input">{nullHint || "null · auto"}</div>
      )}
    </div>
  );
}

// --- feedback / status -------------------------------------------

// Dual-orbit spinner. `state` swaps the live spin for a tick or cross
// when the run terminates; `size` toggles the in-button "sm" variant.
function Spinner({ state = "running", size = "md" }) {
  const cls = [
    "spinner",
    size === "sm" ? "sm" : size === "lg" ? "lg" : "",
    state === "success" ? "done success" : "",
    state === "failure" ? "done failure" : "",
  ].filter(Boolean).join(" ");
  return (
    <div className={cls} role="status" aria-live="polite" aria-label="Working">
      {size !== "sm" && <span className="seed" />}
      {state === "success" && (
        <span className="glyph-overlay"><CheckIcon size={36} /></span>
      )}
      {state === "failure" && (
        <span className="glyph-overlay"><XIcon size={32} /></span>
      )}
    </div>
  );
}

// Bottom-right status pill. `tone` ∈ "neutral" | "success" | "danger"
// (neutral matches the dark-pill default).
function Toast({ tone = "neutral", children }) {
  const cls = `toast ${tone === "success" ? "success" : tone === "danger" ? "danger" : ""}`;
  return <div className={cls}>{children}</div>;
}

// Run scope, phrased one way everywhere: "3 queries × 2 targets · 16 threads".
// Used by the actions-bar ready state and the run overlay subhead.
function formatScope(queries, targets, threads) {
  const q = `${queries} quer${queries === 1 ? "y" : "ies"}`;
  const t = `${targets} target${targets === 1 ? "" : "s"}`;
  return `${q} × ${t} · ${threads} threads`;
}

// Mono mm:ss (or h:mm:ss for >1h runs) formatter — defensive against
// undefined while polling.
function formatElapsed(seconds) {
  const s = Math.max(0, Math.floor(seconds || 0));
  const h = Math.floor(s / 3600);
  const m = Math.floor((s % 3600) / 60);
  const ss = s % 60;
  const pad = (n) => String(n).padStart(2, "0");
  return h > 0 ? `${h}:${pad(m)}:${pad(ss)}` : `${pad(m)}:${pad(ss)}`;
}

// In-bar live indicator used while a run (or save) is mid-flight. The
// `tone` knob lets us reuse this shell post-run to mirror final state.
function ActivityRail({ label, elapsed, tone = "live" }) {
  const cls = `activity-rail ${tone === "live" ? "" : tone} ${tone === "saving" ? "solid" : ""}`;
  return (
    <span className={cls}>
      <span className="live-dot" aria-hidden />
      <span>{label}</span>
      {elapsed != null && <span className="elapsed">· {formatElapsed(elapsed)}</span>}
    </span>
  );
}

function RunOverlay({ open, status, onCancel, onClose, threads, queries, targets }) {
  const state = status?.state || "running";
  const elapsed = status?.elapsed || 0;
  const code = status?.returncode;
  const running = state === "running";
  const logTail = (state === "failed" || state === "cancelled") ? status?.log_tail : null;

  const cardRef = useRef(null);
  const logRef = useRef(null);
  // Esc dismisses only on a terminal state; while running, the only exit is Cancel.
  useDialog(cardRef, open, { onEscape: running ? null : onClose });

  // Snakemake prints the failing rule and traceback at the very end, so jump
  // the log panel to its tail the moment it appears.
  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight;
  }, [logTail]);

  const VIEWS = {
    running: {
      glyph: "running",
      headline: "Hybridization in progress",
      sub: `Running snakemake · ${formatScope(queries, targets, threads)}`,
    },
    succeeded: {
      glyph: "success",
      headline: "Pipeline complete",
      sub: "Snakemake finished cleanly. Outputs are in your results directory.",
    },
    failed: {
      glyph: "failure",
      headline: "Pipeline failed",
      sub: `Snakemake exited with code ${code}. The end of the run log is below.`,
    },
    cancelled: {
      glyph: "failure",
      headline: "Pipeline cancelled",
      sub: "The run was stopped before it finished.",
    },
  };
  const view = VIEWS[state] || VIEWS.running;

  return (
    <div className={`run-overlay ${open ? "open" : ""}`} aria-hidden={!open}>
      <div className="run-card" role="dialog" aria-modal="true" aria-label={view.headline} tabIndex={-1} ref={cardRef}>
        <div className="glyph-wrap">
          <Spinner state={view.glyph} size={running ? "lg" : "md"} />
        </div>
        <div className="ttl">{view.headline}</div>
        <div className="sub">{view.sub}</div>
        <div className="meta">
          <span>elapsed <span className="v">{formatElapsed(elapsed)}</span></span>
          {status?.pid && <span>pid <span className="v">{status.pid}</span></span>}
          {!running && code != null && (
            <span>exit <span className="v">{code}</span></span>
          )}
        </div>
        {logTail && (
          <div className="run-log-wrap">
            <pre className="run-log" ref={logRef} tabIndex={0} aria-label="Pipeline log, last lines">{logTail}</pre>
            {status?.log_path && <div className="run-log-path">Full log: {status.log_path}</div>}
          </div>
        )}
        <div className="ctas">
          <button
            className={running ? "btn sm" : "btn primary sm"}
            onClick={running ? onCancel : onClose}
            aria-label={running ? "Cancel run" : "Dismiss"}
          >
            {running
              ? <React.Fragment><StopIcon /> Cancel run</React.Fragment>
              : <React.Fragment><CheckIcon size={12} /> Dismiss</React.Fragment>}
          </button>
        </div>
      </div>
    </div>
  );
}

// --- layout ------------------------------------------------------

function Hero({ cfg }) {
  const { set, total } = window.countOptionalParams(cfg);
  return (
    <div className="hero">
      <div className="hero-body">
        <h1 className="ttl">RNA → bacterial target hybridization</h1>
        <p className="sub">Build your pipeline config, then hit run.</p>
      </div>
      <div className="hero-meter" title="Optional parameters you've set; the rest stay automatic">
        <span className="v">{set}<span className="sep">/</span>{total}</span>
        <span className="l">optional params set</span>
      </div>
    </div>
  );
}

// Sticky footer holding the readiness summary and primary actions.
// `pipelineState`/`elapsed` come from the parent state machine and
// drive the activity rail when the pipeline is alive.
function ActionsBar({
  cfg, onSave, onRun, onToggleYAML, theme, onTheme,
  pipelineState, elapsed, savingState,
}) {
  const queriesOk = cfg.queries.length >= 1 && cfg.queries.every((q) => q.key.trim() && q.path.trim());
  const targetsOk = cfg.targets.length >= 1 && cfg.targets.every((t) => t.key.trim() && t.path && t.path.trim());
  const ready = queriesOk && targetsOk && cfg.results_dir.trim();

  const running = pipelineState === "running";
  const saving = savingState === "saving";
  const canRun = ready && pipelineState === "idle" && !saving;

  let status;
  if (running) {
    status = <ActivityRail label="Pipeline running" elapsed={elapsed} />;
  } else if (saving) {
    status = <ActivityRail label="Saving config" tone="saving" />;
  } else if (ready) {
    status = (
      <div className="status">
        <span className="status-dot" />
        <span>Ready · {formatScope(cfg.queries.length, cfg.targets.length, cfg.threads)}</span>
      </div>
    );
  } else {
    status = (
      <div className="status">
        <span className="status-dot warn" />
        <span>Resolve required fields to run</span>
      </div>
    );
  }

  return (
    <div className="actions-bar">
      {status}
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
          <span className="accent-dot" />
          <CodeIcon /> config.yaml
        </button>
        <button className="btn" onClick={onSave} data-busy={saving} disabled={running}>
          {saving ? <Spinner size="sm" /> : <SaveIcon />}
          {saving ? "Saving…" : "Save"}
        </button>
        <button
          className="btn primary"
          onClick={onRun}
          disabled={!canRun}
          data-busy={running}
        >
          {running ? <Spinner size="sm" /> : <PlayIcon />}
          {running ? "Running…" : "Run pipeline"}
        </button>
      </div>
    </div>
  );
}

// Slide-up YAML drawer with syntax highlighting and copy/save actions.
function YAMLDrawer({ open, onClose, cfg, onSave }) {
  const lines = useMemo(() => window.buildYAML(cfg), [cfg]);
  const html = useMemo(() => window.renderYAML(lines), [lines]);
  // `plain` is only needed when the user clicks Copy — compute lazily so the
  // drawer doesn't pay for it on every keystroke that mutates `cfg`.
  const copy = () => navigator.clipboard?.writeText(window.renderYAMLPlain(lines));

  const drawerRef = useRef(null);
  useDialog(drawerRef, open, { onEscape: onClose });

  return (
    <React.Fragment>
      <div className={`yaml-backdrop ${open ? "open" : ""}`} onClick={onClose} />
      <aside
        className={`yaml-drawer ${open ? "open" : ""}`}
        aria-hidden={!open}
        role="dialog"
        aria-modal="true"
        aria-label="Generated config"
        tabIndex={-1}
        ref={drawerRef}
      >
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

// export to window for other scripts
Object.assign(window, {
  Toggle, Subcard, PillGroup, PillToggle, SegmentedControl, ErrorBoundary,
  PathInput, KeyedFileRow,
  FolderIcon, PlusIcon, TrashIcon, CodeIcon, PlayIcon, SaveIcon, SunIcon, MoonIcon, CopyIcon,
  CheckIcon, XIcon, StopIcon,
  NullableField, useApiList, useListEditor, promptBrowse,
  Spinner, Toast, ActivityRail, RunOverlay,
  Hero, ActionsBar, YAMLDrawer,
  formatElapsed, formatScope,
});
