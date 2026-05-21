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

// 3. Shared parameters (apply to both RNAhybrid and IntaRNA)
function SharedParamsSection({ cfg, setCfg }) {
  const set = (patch) => setCfg({ ...cfg, ...patch });
  const setNullable = (key, next) => set({ [key]: next });
  return (
    <SectionCard
      num="3"
      title="Shared parameters"
      sub="Apply to both RNAhybrid and IntaRNA"
      anchor="shared"
    >
      <div className="grid-2">
        <Subcard label="Threads">
          <input
            className="input mono"
            type="number"
            min={1}
            step="1"
            value={cfg.threads}
            onChange={(e) => set({ threads: Number(e.target.value) })}
          />
        </Subcard>
        <NullableField name="Max suboptimal hits" kind="integer"
          min={1}
          value={cfg.max_suboptimal_hits}
          onChange={(v) => setNullable("max_suboptimal_hits", v)} />
        <NullableField name="Max total energy" doc="kcal/mol"
          max={0}
          value={cfg.max_total_energy}
          onChange={(v) => setNullable("max_total_energy", v)} />
        <NullableField name="Seed" doc="nucleotides <start>,<end>" kind="pair"
          min={1}
          pairLabels={["start nt", "end nt"]}
          value={cfg.seed}
          onChange={(v) => setNullable("seed", v)} />
      </div>
    </SectionCard>
  );
}

// 4. RNAcalibrate (variant-based)
function RNACalibrateSection({ cfg, setCfg }) {
  const r = cfg.rnacalibrate;
  const set = (patch) => setCfg({ ...cfg, rnacalibrate: { ...r, ...patch } });
  const inert = r.calibration_variant === "uncalibrated";
  return (
    <SectionCard
      num="4"
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
        <Subcard label="Calibration variant">
          <SegmentedControl
            options={window.CALIBRATE_MODES}
            value={r.calibration_variant}
            onChange={(v) => set({ calibration_variant: v })}
            ariaLabel="Calibration variant"
          />
        </Subcard>
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
        <fieldset className="group-fieldset" disabled={inert}>
         <div className="stack lg">
          <div className="grid-2">
            <Subcard label="k (sample size)">
              <input
                className="input mono"
                type="number"
                step="1"
                min={0}
                value={r.k}
                onChange={(e) => set({ k: Number(e.target.value) })}
              />
            </Subcard>
            <Subcard label="Maximum target length" hint="nt">
              <input
                className="input mono"
                type="number"
                step="1"
                min={0}
                value={r.max_target_length}
                onChange={(e) => set({ max_target_length: Number(e.target.value) })}
              />
            </Subcard>
          </div>
          <Subcard
            label="Randomize targets"
            hint="shuffle target sequences before calibration"
            right={
              <Toggle
                on={r.randomize_targets}
                onChange={(v) => set({ randomize_targets: v })}
                ariaLabel="Randomize targets"
              />
            }
          />
         </div>
        </fieldset>
      </div>
    </SectionCard>
  );
}

// 5. RNAhybrid
function RNAHybridSection({ cfg, setCfg }) {
  const r = cfg.rnahybrid;
  const set = (patch) => setCfg({ ...cfg, rnahybrid: { ...r, ...patch } });
  const setF = (key, next) => set({ [key]: next });
  const distLocked = window.isDistributionForced(cfg);

  return (
    <SectionCard
      num="5"
      title="RNAhybrid"
      sub="Target prediction parameters"
      anchor="rnahybrid"
    >
      <div className="stack lg">
        {distLocked ? (
          <div className="field locked">
            <div className="field-label">
              <span>Species model</span>
              <span className="auto-tag locked">LOCKED · overridden</span>
            </div>
            <div className="ghost-input">
              overridden · distribution is provided by RNAcalibrate
            </div>
          </div>
        ) : (
          <Subcard label="Species model" hint="3′UTR distribution">
            <select
              className="select"
              value={r.species}
              onChange={(e) => set({ species: e.target.value })}
            >
              {window.SPECIES_OPTIONS.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
        )}

        <div className="grid-2">
          <NullableField name="Max internal loop" kind="integer"
            min={0}
            value={r.max_internal_loop} onChange={(v) => setF("max_internal_loop", v)} />
          <NullableField name="Max bulge loop" kind="integer"
            min={0}
            value={r.max_bulge_loop} onChange={(v) => setF("max_bulge_loop", v)} />
          <NullableField name="p-value threshold"
            min={0}
            value={r.pvalue_threshold} onChange={(v) => setF("pvalue_threshold", v)} />
          <NullableField name="Distribution" doc="<mean>,<std>" kind="pair"
            min={0}
            pairLabels={["mean", "std"]}
            value={r.distribution}
            onChange={(v) => setF("distribution", v)}
            forcedNull={distLocked}
            forcedNullHint="null · distribution is provided by RNAcalibrate"
          />
        </div>
      </div>
    </SectionCard>
  );
}

// 6. IntaRNA
function IntaRNASection({ cfg, setCfg }) {
  const it = cfg.intarna;
  const set     = (patch) => setCfg({ ...cfg, intarna: { ...it, ...patch } });
  const setSeed = (patch) => set({ seed:          { ...it.seed,          ...patch } });
  const setAcc  = (patch) => set({ accessibility: { ...it.accessibility, ...patch } });
  const setOut  = (patch) => set({ output:        { ...it.output,        ...patch } });

  const accInert     = it.accessibility_variant === "off";
  const seedDisabled = !it.seed.enabled;
  const seedDerived  = window.isIntarnaSeedDerived(cfg);

  return (
    <SectionCard
      num="6"
      title="IntaRNA"
      sub="Accessibility-aware RNA–RNA interaction prediction"
      anchor="intarna"
    >
      <div className="stack lg">
        <Subcard label="Accessibility variant">
          <SegmentedControl
            options={window.INTARNA_ACC_MODES}
            value={it.accessibility_variant}
            onChange={(v) => set({ accessibility_variant: v })}
            ariaLabel="Accessibility variant"
          />
        </Subcard>

        <div className="grid-2">
          <Subcard label="Prediction mode">
            <select
              className="select"
              value={it.prediction_mode}
              onChange={(e) => set({ prediction_mode: e.target.value })}
            >
              {window.INTARNA_PREDICTION_MODES.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
          <Subcard label="Interaction model">
            <select
              className="select"
              value={it.model}
              onChange={(e) => set({ model: e.target.value })}
            >
              {window.INTARNA_MODELS.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
          <Subcard label="Max interaction length" hint="0 = auto">
            <input
              className="input mono"
              type="number"
              min={0}
              step="1"
              value={it.max_interaction_length}
              onChange={(e) => set({ max_interaction_length: Number(e.target.value) })}
            />
          </Subcard>
          <Subcard label="Max loop size">
            <input
              className="input mono"
              type="number"
              min={0}
              step="1"
              value={it.max_loop_size}
              onChange={(e) => set({ max_loop_size: Number(e.target.value) })}
            />
          </Subcard>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 10 }}>
            <div className="eyebrow strong">Seed</div>
            <div style={{ flex: 1 }} />
            <span style={{ fontSize: 12, color: seedDisabled ? "var(--fg-4)" : "var(--accent-ink)", marginRight: 8 }}>
              {seedDisabled ? "Disabled" : "Enabled"}
            </span>
            <Toggle on={it.seed.enabled} onChange={(v) => setSeed({ enabled: v })} ariaLabel="Enable seed enforcement" />
          </div>
          <fieldset className="group-fieldset" disabled={seedDisabled}>
           <div className="stack sm">
            <div className="grid-2">
              <NullableField name="Length" doc="[2..20]" kind="integer"
                min={2} max={20}
                value={it.seed.length} onChange={(v) => setSeed({ length: v })}
                forcedNull={seedDerived}
                forcedNullHint="null · derived from shared seed" />
              <NullableField name="Query range" doc="'from-to,...'"
                value={it.seed.query_range} onChange={(v) => setSeed({ query_range: v })}
                forcedNull={seedDerived}
                forcedNullHint="null · derived from shared seed" />
              <NullableField name="Max unpaired bases" kind="integer"
                min={0}
                value={it.seed.max_unpaired_bases} onChange={(v) => setSeed({ max_unpaired_bases: v })}
                forcedNull={seedDerived}
                forcedNullHint="null · derived from shared seed (set to 0)" />
              <NullableField name="Target range" doc="'from-to,...'"
                value={it.seed.target_range} onChange={(v) => setSeed({ target_range: v })} />
              <NullableField name="Max energy"
                max={999}
                value={it.seed.max_energy} onChange={(v) => setSeed({ max_energy: v })} />
              <NullableField name="Max hybrid energy"
                max={999}
                value={it.seed.max_hybrid_energy} onChange={(v) => setSeed({ max_hybrid_energy: v })} />
              <NullableField name="Min unpaired probability" doc="[0..1]"
                min={0} max={1}
                value={it.seed.min_unpaired_probability} onChange={(v) => setSeed({ min_unpaired_probability: v })} />
            </div>
            <PillGroup label="Constraints">
              <PillToggle on={it.seed.forbid_gu}         onChange={(v) => setSeed({ forbid_gu: v })}         label="Forbid G:U pairs" />
              <PillToggle on={it.seed.forbid_gu_at_ends} onChange={(v) => setSeed({ forbid_gu_at_ends: v })} label="Forbid G:U at ends" />
              <PillToggle on={it.seed.report_best_only}  onChange={(v) => setSeed({ report_best_only: v })}  label="Report best only" />
            </PillGroup>
           </div>
          </fieldset>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 10 }}>
            <div className="eyebrow strong">Accessibility</div>
            {accInert && (
              <span style={{ fontSize: 11, color: "var(--fg-4)", marginLeft: 8, fontStyle: "italic" }}>
                · variant is Off — flags below are sent but ignored
              </span>
            )}
          </div>
          <fieldset className="group-fieldset" disabled={accInert}>
           <div className="stack sm">
            <div className="grid-2">
              <NullableField name="Window" doc="0 = global, null = default" kind="integer"
                min={0}
                value={it.accessibility.window} onChange={(v) => setAcc({ window: v })} />
              <NullableField name="Max bp span" doc="0 = global, null = default" kind="integer"
                min={0}
                value={it.accessibility.max_bp_span} onChange={(v) => setAcc({ max_bp_span: v })} />
            </div>
            <PillGroup label="Constraints">
              <PillToggle on={it.accessibility.forbid_lonely_pairs} onChange={(v) => setAcc({ forbid_lonely_pairs: v })} label="Forbid lonely pairs" />
              <PillToggle on={it.accessibility.forbid_gu_at_ends}   onChange={(v) => setAcc({ forbid_gu_at_ends: v })}   label="Forbid G:U at ends" />
            </PillGroup>
           </div>
          </fieldset>
        </div>

        <div>
          <div className="eyebrow strong" style={{ marginBottom: 10 }}>Output</div>
          <div className="grid-2">
            <NullableField name="Max ΔE" doc="cap vs MFE"
              value={it.output.max_delta_energy} onChange={(v) => setOut({ max_delta_energy: v })} />
            <Subcard label="Overlap">
              <select
                className="select"
                value={it.output.overlap}
                onChange={(e) => setOut({ overlap: e.target.value })}
              >
                {window.INTARNA_OVERLAPS.map((o) => (
                  <option key={o.value} value={o.value}>{o.label}</option>
                ))}
              </select>
            </Subcard>
            <NullableField name="Min unpaired probability" doc="[0..1]"
              min={0} max={1}
              value={it.output.min_unpaired_probability} onChange={(v) => setOut({ min_unpaired_probability: v })} />
            <PillGroup label="Constraints">
              <PillToggle on={it.output.forbid_lonely_pairs} onChange={(v) => setOut({ forbid_lonely_pairs: v })} label="Forbid lonely pairs" />
              <PillToggle on={it.output.forbid_gu_at_ends}   onChange={(v) => setOut({ forbid_gu_at_ends: v })}   label="Forbid G:U at ends" />
            </PillGroup>
          </div>
        </div>
      </div>
    </SectionCard>
  );
}

// 7. Plots
function PlotsSection({ cfg, setCfg }) {
  const pl = cfg.plots;
  const set = (patch) => setCfg({ ...cfg, plots: { ...pl, ...patch } });
  const enabled = pl.types.length > 0;
  const setEnabled = (v) => set({ types: v ? window.PLOT_TYPES.map((p) => p.value) : [] });
  const toggleType = (val) => {
    set({ types: pl.types.includes(val) ? pl.types.filter((t) => t !== val) : [...pl.types, val] });
  };

  const loci = useListEditor(
    pl.locus,
    (next) => set({ locus: next }),
    () => ({ id: "l" + Date.now().toString(36), value: "", mirna: "" }),
  );
  const genes = useListEditor(
    pl.gene,
    (next) => set({ gene: next }),
    () => ({ id: "g" + Date.now().toString(36), value: "", mirna: "" }),
  );
  const proteins = useListEditor(
    pl.protein,
    (next) => set({ protein: next }),
    () => ({ id: "p" + Date.now().toString(36), value: "", mirna: "" }),
  );

  return (
    <SectionCard
      num="7"
      title="Plots"
      sub="Optional · omit or disable to skip"
      anchor="plots"
      right={
        <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
          <span style={{ fontSize: 12, color: enabled ? "var(--accent-ink)" : "var(--fg-4)" }}>
            {enabled ? "Enabled" : "Disabled"}
          </span>
          <Toggle on={enabled} onChange={setEnabled} ariaLabel="Enable plots" />
        </div>
      }
    >
      <fieldset className="group-fieldset" disabled={!enabled}>
       <div className="stack lg">
        <Subcard label="Plot types" hint={`${pl.types.length}/${window.PLOT_TYPES.length} on`}>
          <div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>
            {window.PLOT_TYPES.map((p) => {
              const on = pl.types.includes(p.value);
              return (
                <PillToggle key={p.value} on={on} onChange={() => toggleType(p.value)} label={p.label} />
              );
            })}
          </div>
        </Subcard>
        <div className="grid-3">
          <Subcard label="Base font size">
            <input
              className="input mono"
              type="number"
              step="1"
              min={0}
              value={pl.basesize}
              onChange={(e) => set({ basesize: Number(e.target.value) })}
            />
          </Subcard>
          <Subcard label="p-value threshold" hint="0 – 1">
            <input
              className="input mono"
              type="number"
              step="0.001"
              min={0}
              value={pl.pvalue_threshold}
              onChange={(e) => set({ pvalue_threshold: Number(e.target.value) })}
            />
          </Subcard>
          <Subcard label="Top-N miRNAs" hint="per-miRNA plot cap">
            <input
              className="input mono"
              type="number"
              min={1}
              step="1"
              value={pl.per_mirna_top_n}
              onChange={(e) => set({ per_mirna_top_n: Number(e.target.value) })}
            />
          </Subcard>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 8 }}>
            <div className="eyebrow strong">Locus tags · {pl.locus.length}</div>
            <div style={{ flex: 1 }} />
            <button className="btn sm" onClick={loci.add}><PlusIcon /> Add locus</button>
          </div>
          <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8, padding: "0 2px", marginBottom: 6 }}>
            <span>Locus tag</span>
            <span>miRNA <span style={{ textTransform: "none", fontWeight: 400 }}>· optional</span></span>
            <span></span>
          </div>
          <div className="stack xs">
            {pl.locus.map((l, i) => {
              const orphanMirna = !!(l.mirna && l.mirna.trim() && !(l.value && l.value.trim()));
              return (
                <div key={l.id} style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8 }}>
                  <input
                    className="input mono"
                    value={l.value}
                    placeholder="locus tag (required), e.g. FE838_RS16060"
                    onChange={(e) => loci.update(i, { value: e.target.value })}
                    style={{ borderColor: orphanMirna ? "var(--danger)" : undefined }}
                    aria-invalid={orphanMirna}
                  />
                  <input
                    className="input mono"
                    value={l.mirna}
                    placeholder="miRNA (optional), e.g. hsa-miR-1226-5p"
                    onChange={(e) => loci.update(i, { mirna: e.target.value })}
                  />
                  <button className="btn sm danger-ghost" onClick={() => loci.remove(i)} aria-label="Remove locus">
                    <TrashIcon />
                  </button>
                </div>
              );
            })}
            {pl.locus.length === 0 && (
              <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic" }}>No locus tags — list will be empty.</div>
            )}
          </div>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 8 }}>
            <div className="eyebrow strong">Genes · {pl.gene.length}</div>
            <div style={{ flex: 1 }} />
            <button className="btn sm" onClick={genes.add}><PlusIcon /> Add gene</button>
          </div>
          <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8, padding: "0 2px", marginBottom: 6 }}>
            <span>Gene</span>
            <span>miRNA <span style={{ textTransform: "none", fontWeight: 400 }}>· optional</span></span>
            <span></span>
          </div>
          <div className="stack xs">
            {pl.gene.map((g, i) => {
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
            {pl.gene.length === 0 && (
              <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic" }}>No genes — list will be empty.</div>
            )}
          </div>
        </div>

        <div>
          <div style={{ display: "flex", alignItems: "center", marginBottom: 8 }}>
            <div className="eyebrow strong">Proteins · {pl.protein.length}</div>
            <div style={{ flex: 1 }} />
            <button className="btn sm" onClick={proteins.add}><PlusIcon /> Add protein</button>
          </div>
          <div className="eyebrow" style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8, padding: "0 2px", marginBottom: 6 }}>
            <span>Protein</span>
            <span>miRNA <span style={{ textTransform: "none", fontWeight: 400 }}>· optional</span></span>
            <span></span>
          </div>
          <div className="stack xs">
            {pl.protein.map((p, i) => {
              const orphanMirna = !!(p.mirna && p.mirna.trim() && !(p.value && p.value.trim()));
              return (
                <div key={p.id} style={{ display: "grid", gridTemplateColumns: "180px 1fr auto", gap: 8 }}>
                  <input
                    className="input mono"
                    value={p.value}
                    placeholder="protein name (required)"
                    onChange={(e) => proteins.update(i, { value: e.target.value })}
                    style={{ borderColor: orphanMirna ? "var(--danger)" : undefined }}
                    aria-invalid={orphanMirna}
                  />
                  <input
                    className="input mono"
                    value={p.mirna}
                    placeholder="miRNA (optional), e.g. hsa-miR-1226-5p"
                    onChange={(e) => proteins.update(i, { mirna: e.target.value })}
                  />
                  <button className="btn sm danger-ghost" onClick={() => proteins.remove(i)} aria-label="Remove protein">
                    <TrashIcon />
                  </button>
                </div>
              );
            })}
            {pl.protein.length === 0 && (
              <div style={{ fontSize: 12, color: "var(--fg-4)", fontStyle: "italic" }}>No proteins — list will be empty.</div>
            )}
          </div>
        </div>
       </div>
      </fieldset>
    </SectionCard>
  );
}

// 8. Output
function OutputSection({ cfg, setCfg }) {
  return (
    <SectionCard num="8" title="Output" sub="Destination folder" anchor="output">
      <Subcard label="Results directory" hint="created if missing">
        <PathInput
          value={cfg.results_dir}
          onChange={(v) => setCfg({ ...cfg, results_dir: v })}
          placeholder="Data/Results/"
        />
      </Subcard>
    </SectionCard>
  );
}

Object.assign(window, {
  SectionCard,
  QueriesSection, TargetsSection, SharedParamsSection,
  RNACalibrateSection, RNAHybridSection, IntaRNASection,
  PlotsSection, OutputSection
});
