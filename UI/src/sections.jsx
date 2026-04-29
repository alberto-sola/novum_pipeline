// Novum Pipeline — section UIs

function SectionCard({ num, title, sub, anchor, children, right }) {
  return (
    <section className="card" data-anchor={anchor}>
      <div className="card-head">
        <div className="num">{num}</div>
        <div className="title">{title}</div>
        <div className="sub">{sub}</div>
        <div style={{ flex: 1 }} />
        {right}
      </div>
      <div className="card-body">{children}</div>
    </section>
  );
}

// 1. Queries (miRNA FASTAs, keyed)
function QueriesSection({ cfg, setCfg }) {
  const { items, hasApi } = useApiList("list_queries");
  const { add, update, remove } = useListEditor(
    cfg.queries,
    (next) => setCfg({ ...cfg, queries: next }),
    () => ({ id: "q" + Date.now().toString(36), key: "", path: "" }),
    1,
  );
  const cols = hasApi ? "180px 1fr auto" : "180px 1fr auto auto";
  return (
    <SectionCard
      num="1"
      title="Queries"
      sub="One query FASTA per sample (joined by key)"
      anchor="queries"
      right={<button className="btn sm" onClick={add}><PlusIcon /> Add query</button>}
    >
      <div className="stack sm">
        <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: cols, gap: 8, padding: "0 2px" }}>
          <span>Sample key</span>
          <span>{hasApi ? "Query FASTA · Data/Raw/RNAs/" : "Query FASTA path"}</span>
          <span></span>
          {!hasApi && <span></span>}
        </div>
        {cfg.queries.length === 0 && (
          <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic", padding: "8px 2px" }}>
            No queries yet — click <b>Add query</b>.
          </div>
        )}
        {cfg.queries.map((q, i) => (
          <KeyedFileRow
            key={q.id}
            keyVal={q.key}
            pathVal={q.path}
            items={hasApi ? items : undefined}
            onKey={(v) => update(i, { key: v })}
            onPath={(v) => update(i, { path: v })}
            onRemove={() => remove(i)}
            canRemove={cfg.queries.length > 1}
            placeholderKey="e.g. escherichia"
            placeholderPath="Data/Raw/RNAs/validated_*.fa"
            emptyHint="No files in Data/Raw/RNAs/"
          />
        ))}
        <div style={{ fontSize: 12, color: "var(--fg-4)", paddingTop: 4 }}>
          {cfg.queries.length} quer{cfg.queries.length === 1 ? "y" : "ies"} · keys must match target keys
        </div>
      </div>
    </SectionCard>
  );
}

// 2. Targets — keyed (CDS / RNAs / full genomes)
function TargetsSection({ cfg, setCfg }) {
  const { items, hasApi } = useApiList("list_targets");
  const { add, update, remove } = useListEditor(
    cfg.targets,
    (next) => setCfg({ ...cfg, targets: next }),
    () => ({ id: "t" + Date.now().toString(36), key: "", path: "" }),
    1,
  );
  const queryKeys = cfg.queries.map((q) => q.key).filter(Boolean);
  const cols = hasApi ? "180px 1fr auto" : "180px 1fr auto auto";
  return (
    <SectionCard
      num="2"
      title="Targets"
      sub="One or more (CDS, RNAs, full genomes)"
      anchor="targets"
      right={<button className="btn sm" onClick={add}><PlusIcon /> Add target</button>}
    >
      <div className="stack sm">
        <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: cols, gap: 8, padding: "0 2px" }}>
          <span>Sample key</span>
          <span>{hasApi ? "Target genome · Data/Raw/genomes/" : "Target genome path"}</span>
          <span></span>
          {!hasApi && <span></span>}
        </div>
        {cfg.targets.length === 0 && (
          <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic", padding: "8px 2px" }}>
            No targets yet — click <b>Add target</b>.
          </div>
        )}
        {cfg.targets.map((t, i) => {
          const keyMatched = queryKeys.includes(t.key);
          return (
            <div key={t.id} className="stack xs">
              <KeyedFileRow
                keyVal={t.key}
                pathVal={t.path}
                items={hasApi ? items : undefined}
                onKey={(v) => update(i, { key: v })}
                onPath={(v) => update(i, { path: v })}
                onRemove={() => remove(i)}
                canRemove={cfg.targets.length > 1}
                placeholderKey="e.g. escherichia"
                placeholderPath="Data/Raw/genomes/*.fna"
                emptyHint="No files in Data/Raw/genomes/"
              />
              <div style={{ fontSize: 12, color: t.key && !keyMatched ? "var(--warn)" : "var(--fg-4)", padding: "2px 4px" }}>
                {t.key
                  ? (keyMatched
                      ? <span>✓ matches query <code style={{ fontFamily: "var(--font-mono)" }}>{t.key}</code></span>
                      : <span>No query named <code style={{ fontFamily: "var(--font-mono)" }}>{t.key}</code></span>)
                  : <span>Enter a key shared with a query</span>}
              </div>
            </div>
          );
        })}
        <div style={{ fontSize: 12, color: "var(--fg-4)", paddingTop: 4 }}>
          {cfg.targets.length} target{cfg.targets.length === 1 ? "" : "s"}
        </div>
      </div>
    </SectionCard>
  );
}

// 3. RNAcalibrate (mode-based)
function RNACalibrateSection({ cfg, setCfg }) {
  const r = cfg.rnacalibrate;
  const set = (patch) => setCfg({ ...cfg, rnacalibrate: { ...r, ...patch } });
  const inert = r.mode === "uncalibrated";
  return (
    <SectionCard
      num="3"
      title="RNAcalibrate"
      sub="Calibrate hybridization against random targets"
      anchor="rnacalibrate"
      right={
        inert ? (
          <span style={{
            fontSize: 11, letterSpacing: "0.04em", textTransform: "uppercase",
            fontWeight: 600, color: "var(--fg-4)",
            padding: "3px 8px", border: "1px solid var(--line)", borderRadius: 999,
            background: "var(--surface-2)"
          }}>
            Parameters inactive
          </span>
        ) : null
      }
    >
      <div className="stack lg">
        <Field label="Mode" hint="Calibrated · Uncalibrated · Both">
          <div className="variation-switch" role="tablist" style={{ display: "inline-flex" }}>
            {window.CALIBRATE_MODES.map((m) => (
              <button
                key={m.value}
                className={r.mode === m.value ? "active" : ""}
                onClick={() => set({ mode: m.value })}
                style={{ padding: "6px 14px" }}
              >
                {m.label}
              </button>
            ))}
          </div>
        </Field>
        {inert && (
          <div style={{
            fontSize: 12.5, color: "var(--fg-3)",
            padding: "10px 12px",
            background: "var(--surface-2)",
            border: "1px dashed var(--line)",
            borderRadius: "var(--r-sm)",
            display: "flex", gap: 10, alignItems: "flex-start"
          }}>
            <span aria-hidden style={{
              flex: "0 0 auto", width: 18, height: 18, borderRadius: "50%",
              background: "var(--surface)", border: "1px solid var(--line)",
              display: "inline-flex", alignItems: "center", justifyContent: "center",
              fontSize: 11, fontWeight: 700, color: "var(--fg-3)"
            }}>i</span>
            <div>
              Calibration is off. The parameters below won’t affect the run — set
              RNAhybrid’s <code style={{ fontFamily: "var(--font-mono)" }}>distribution</code> instead.
            </div>
          </div>
        )}
        <div
          className="stack lg"
          style={{
            opacity: inert ? 0.42 : 1,
            filter: inert ? "saturate(0.2)" : "none",
            transition: "opacity 0.15s, filter 0.15s",
            pointerEvents: inert ? "none" : "auto"
          }}
          aria-disabled={inert}
        >
          <div className="grid-2">
            <Field label="k (sample size)" hint="int">
              <input
                className="input mono"
                type="number"
                step="1"
                value={r.k}
                disabled={inert}
                onChange={(e) => set({ k: Number(e.target.value) })}
              />
            </Field>
            <Field label="Maximum target length" hint="int · nt">
              <input
                className="input mono"
                type="number"
                step="1"
                value={r.max_target_length}
                disabled={inert}
                onChange={(e) => set({ max_target_length: Number(e.target.value) })}
              />
            </Field>
          </div>
          <div style={{
            display: "flex", alignItems: "center", justifyContent: "space-between",
            padding: "10px 12px", background: "var(--surface-2)", border: "1px solid var(--line)",
            borderRadius: "var(--r-sm)"
          }}>
            <div>
              <div style={{ fontSize: 13, fontWeight: 500 }}>Randomize targets</div>
              <div style={{ fontSize: 12, color: "var(--fg-3)", marginTop: 2 }}>
                Shuffle target sequences before calibration
              </div>
            </div>
            <Check
              on={r.randomize_targets}
              onChange={(v) => !inert && set({ randomize_targets: v })}
              label=""
            />
          </div>
        </div>
      </div>
    </SectionCard>
  );
}

// 4. RNAhybrid
function RNAHybridSection({ cfg, setCfg }) {
  const r = cfg.rnahybrid;
  const set = (patch) => setCfg({ ...cfg, rnahybrid: { ...r, ...patch } });
  const setF = (key, next) => set({ [key]: next });
  const distLocked = window.isDistributionForced(cfg);

  return (
    <SectionCard
      num="4"
      title="RNAhybrid"
      sub="Target prediction parameters"
      anchor="rnahybrid"
    >
      <div className="stack lg">
        <div className="grid-2">
          <Field label="Threads" hint="int">
            <input
              className="input mono"
              type="number"
              min={1}
              step="1"
              value={r.threads}
              onChange={(e) => set({ threads: Number(e.target.value) })}
            />
          </Field>
          <div className={`field ${distLocked ? "locked" : ""}`}>
            <div className="field-label">
              <span>Species model</span>
              {distLocked ? (
                <span className="auto-tag locked">LOCKED · overridden</span>
              ) : (
                <span className="hint">3′UTR distribution</span>
              )}
            </div>
            {distLocked ? (
              <div className="ghost-input">
                overridden · distribution is provided by RNAcalibrate
              </div>
            ) : (
              <select
                className="select"
                value={r.species}
                onChange={(e) => set({ species: e.target.value })}
              >
                {window.SPECIES_OPTIONS.map((o) => (
                  <option key={o.value} value={o.value}>{o.label}</option>
                ))}
              </select>
            )}
          </div>
        </div>

        <div>
          <div className="eyebrow strong" style={{ marginBottom: 10 }}>
            Optional parameters — toggle to set, else null
          </div>
          <div className="grid-2">
            <NullableField name="Energy" doc="min free energy (kcal/mol)"
              value={r.energy} onChange={(v) => setF("energy", v)} />
            <NullableField name="p-value" doc="p-value cutoff"
              value={r.pvalue} onChange={(v) => setF("pvalue", v)} />
            <NullableField name="Hits" doc="max hits per target" kind="integer"
              value={r.hits} onChange={(v) => setF("hits", v)} />
            <NullableField name="U" doc="max internal loop size" kind="integer"
              value={r.u} onChange={(v) => setF("u", v)} />
            <NullableField name="V" doc="max bulge loop size" kind="integer"
              value={r.v} onChange={(v) => setF("v", v)} />
            <NullableField name="Seed" doc="start,end nt" kind="pair"
              pairLabels={["start nt", "end nt"]}
              value={r.seed} onChange={(v) => setF("seed", v)} />
            <div style={{ gridColumn: "1 / -1" }}>
              <NullableField name="distribution" doc="mean,std (integers)" kind="pair"
                pairLabels={["mean", "std"]}
                value={r.distribution}
                onChange={(v) => setF("distribution", v)}
                forcedNull={distLocked}
                forcedNullHint="null · distribution is provided by RNAcalibrate"
              />
            </div>
          </div>
        </div>
      </div>
    </SectionCard>
  );
}

// 5. Build plots
function BuildPlotsSection({ cfg, setCfg }) {
  const bp = cfg.build_plots;
  const set = (patch) => setCfg({ ...cfg, build_plots: { ...bp, ...patch } });
  const allPlotTypes = window.PLOT_TYPES.map((p) => p.value);
  const enabled = bp.types.length > 0;
  const setEnabled = (v) => set({ types: v ? allPlotTypes : [] });
  const toggleType = (val) => {
    set({ types: bp.types.includes(val) ? bp.types.filter((t) => t !== val) : [...bp.types, val] });
  };

  const loci = useListEditor(
    bp.locus,
    (next) => set({ locus: next }),
    () => ({ id: "l" + Date.now().toString(36), value: "" }),
  );
  const genes = useListEditor(
    bp.gene,
    (next) => set({ gene: next }),
    () => ({ id: "g" + Date.now().toString(36), value: "", mirna: "" }),
  );

  return (
    <SectionCard
      num="5"
      title="Build plots"
      sub="Optional · omit or disable to skip"
      anchor="build_plots"
      right={
        <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
          <span style={{ fontSize: 12, color: enabled ? "var(--accent-ink)" : "var(--fg-4)" }}>
            {enabled ? "Enabled" : "Disabled"}
          </span>
          <Toggle on={enabled} onChange={setEnabled} ariaLabel="Enable build_plots" />
        </div>
      }
    >
      <div className="stack lg" style={{ opacity: enabled ? 1 : 0.5, pointerEvents: enabled ? "auto" : "none", transition: "opacity 0.15s" }}>
        <Field label="Plot types" hint="select any · all three = “all”">
          <div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>
            {window.PLOT_TYPES.map((p) => {
              const on = bp.types.includes(p.value);
              return (
                <button
                  key={p.value}
                  type="button"
                  onClick={() => toggleType(p.value)}
                  className="btn sm"
                  style={{
                    background: on ? "var(--accent-soft)" : "var(--surface)",
                    color: on ? "var(--accent-ink)" : "var(--fg-2)",
                    borderColor: on ? "var(--accent)" : "var(--line)"
                  }}
                >
                  {on ? "✓ " : ""}{p.label}
                </button>
              );
            })}
          </div>
        </Field>
        <div className="grid-2">
          <Field label="Base font size" hint="int">
            <input
              className="input mono"
              type="number"
              step="1"
              value={bp.basesize}
              onChange={(e) => set({ basesize: Number(e.target.value) })}
            />
          </Field>
          <Field label="p-value cutoff" hint="float · 0 – 1">
            <input
              className="input mono"
              type="number"
              step="0.001"
              value={bp.pvalue}
              onChange={(e) => set({ pvalue: Number(e.target.value) })}
            />
          </Field>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 8 }}>
            <div className="eyebrow strong">Locus tags · {bp.locus.length}</div>
            <div style={{ flex: 1 }} />
            <button className="btn sm" onClick={loci.add}><PlusIcon /> Add locus</button>
          </div>
          <div className="stack xs">
            {bp.locus.map((l, i) => (
              <div key={l.id} style={{ display: "grid", gridTemplateColumns: "1fr auto", gap: 8 }}>
                <input
                  className="input mono"
                  value={l.value}
                  placeholder="e.g. FE838_RS16060"
                  onChange={(e) => loci.update(i, { value: e.target.value })}
                />
                <button className="btn sm danger-ghost" onClick={() => loci.remove(i)} aria-label="Remove locus">
                  <TrashIcon />
                </button>
              </div>
            ))}
            {bp.locus.length === 0 && (
              <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic" }}>No locus tags — list will be empty.</div>
            )}
          </div>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 8 }}>
            <div className="eyebrow strong">Genes · {bp.gene.length}</div>
            <div style={{ flex: 1 }} />
            <button className="btn sm" onClick={genes.add}><PlusIcon /> Add gene</button>
          </div>
          <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8, padding: "0 2px", marginBottom: 6 }}>
            <span>Gene</span>
            <span>miRNA <span style={{ textTransform: "none", fontWeight: 400 }}>· optional</span></span>
            <span></span>
          </div>
          <div className="stack xs">
            {bp.gene.map((g, i) => {
              const orphanMirna = !!(g.mirna && g.mirna.trim() && !(g.value && g.value.trim()));
              return (
                <div key={g.id} style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8 }}>
                  <input
                    className="input mono"
                    value={g.value}
                    placeholder="gene name (required), e.g. yegH"
                    onChange={(e) => genes.update(i, { value: e.target.value })}
                    style={{ borderColor: orphanMirna ? "var(--danger)" : undefined }}
                    aria-invalid={orphanMirna}
                  />
                  <input
                    className="input mono"
                    value={g.mirna}
                    placeholder="miRNA (optional), e.g. hsa-miR-1226-5p"
                    onChange={(e) => genes.update(i, { mirna: e.target.value })}
                  />
                  <button className="btn sm danger-ghost" onClick={() => genes.remove(i)} aria-label="Remove gene">
                    <TrashIcon />
                  </button>
                </div>
              );
            })}
            {bp.gene.length === 0 && (
              <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic" }}>No genes — list will be empty.</div>
            )}
          </div>
        </div>
      </div>
    </SectionCard>
  );
}

// 6. Output
function OutputSection({ cfg, setCfg }) {
  return (
    <SectionCard num="6" title="Output" sub="Destination folder" anchor="output">
      <Field label="Results directory" hint="created if missing">
        <PathInput
          value={cfg.results_dir}
          onChange={(v) => setCfg({ ...cfg, results_dir: v })}
          placeholder="Data/Results/"
        />
      </Field>
    </SectionCard>
  );
}

Object.assign(window, {
  SectionCard,
  QueriesSection, TargetsSection, RNACalibrateSection, RNAHybridSection, BuildPlotsSection, OutputSection
});
