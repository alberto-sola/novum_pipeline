// Novum Pipeline — shared form components (React, Babel)

const { useState, useRef, useEffect, useMemo } = React;

// --- primitives ---------------------------------------------------

function Toggle({ on, onChange, ariaLabel }) {
  return (
    <button
      type="button"
      className={`toggle ${on ? "on" : ""}`}
      onClick={() => onChange(!on)}
      aria-label={ariaLabel}
      aria-pressed={on}
    />
  );
}

function Check({ on, onChange, label }) {
  return (
    <label className={`check ${on ? "on" : ""}`}>
      <input type="checkbox" className="sr-only" checked={on} onChange={(e) => onChange(e.target.checked)} />
      <span className="check-box" />
      <span>{label}</span>
    </label>
  );
}

function Field({ label, hint, children }) {
  return (
    <div className="field">
      <div className="field-label">
        <span>{label}</span>
        {hint && <span className="hint">{hint}</span>}
      </div>
      {children}
    </div>
  );
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

// --- nullable numeric input --------------------------------------

function NullableField({
  name, doc, value, onChange,
  kind = "number",        // "number" | "integer" | "pair"
  placeholder,
  pairLabels,             // ["start nt", "end nt"] etc.
  forcedNull = false,
  forcedNullHint,
  min,
  max
}) {
  const on = !forcedNull && value.set;
  const toInt = (s) => {
    if (s === "" || s === "-") return s;
    const n = parseInt(s, 10);
    return Number.isFinite(n) ? n : "";
  };
  return (
    <div className={`nullable ${on ? "on" : ""} ${forcedNull ? "locked" : ""}`}>
      <div className="top">
        <div className="name">
          {name} <em>· {doc}</em>
        </div>
        <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
          <span className="auto-tag">
            {forcedNull ? "LOCKED · null" : (on ? "SET" : "AUTO (null)")}
          </span>
          <Toggle
            on={on}
            onChange={(v) => !forcedNull && onChange({ ...value, set: v })}
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
              const v = kind === "integer" ? toInt(raw) : Number(raw);
              onChange({ ...value, value: v });
            }}
          />
        )
      ) : (
        <div className="ghost-input">null · let RNAhybrid decide</div>
      )}
    </div>
  );
}

// export to window for other scripts
Object.assign(window, {
  Toggle, Check, Field, PathInput, KeyedFileRow,
  FolderIcon, PlusIcon, TrashIcon, CodeIcon, PlayIcon, SaveIcon, SunIcon, MoonIcon, CopyIcon,
  NullableField, useApiList, useListEditor, promptBrowse,
});
