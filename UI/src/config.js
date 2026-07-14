// Novum Pipeline — config <-> form-state bridge.
//
// `configToObject` is the sole writer of the config schema: api.save_config dumps
// it, and yaml.js renders the preview from it.

// Nullable form field -> serialized {v, kind}. Pair-shaped fields (`a`/`b`) emit
// "a,b"; unset/forced/blank fields collapse to null. The one place that collapsing
// happens — data.js's derived warnings must route through it too, or they can
// contradict what Save writes.
function resolveNullable(field, forced) {
  if (forced || !field.set) {
    return { v: null, kind: "null" };
  }
  if ("a" in field && "b" in field) {
    const { a, b } = field;
    if (a === "" || b === "" || a == null || b == null) {
      return { v: null, kind: "null" };
    }
    return { v: `${a},${b}`, kind: "str" };
  }
  if (field.value === "" || field.value === null || field.value === undefined) {
    return { v: null, kind: "null" };
  }
  return { v: field.value, kind: typeof field.value === "string" ? "str" : "num" };
}
window.resolveNullable = resolveNullable;

function _nullableValue(field, forced) {
  return resolveNullable(field, forced).v;
}

function _pairListToArray(items) {
  return (items || [])
    .map((it) => ({ v: (it.value || "").trim(), m: (it.mirna || "").trim() }))
    .filter((it) => it.v)
    .map((it) => (it.m ? `${it.v},${it.m}` : it.v));
}

// form-state cfg -> plain config object (what launcher.save_config/run_pipeline receive)
window.configToObject = function configToObject(cfg) {
  const obj = {};

  obj.queries = {};
  cfg.queries.forEach((q) => { obj.queries[(q.key || "").trim() || "unnamed_query"] = q.path || ""; });
  obj.targets = {};
  cfg.targets.forEach((t) => { obj.targets[(t.key || "").trim() || "unnamed_target"] = t.path || ""; });

  obj.threads = cfg.threads;
  window.SHARED_NULLABLE_KEYS.forEach((k) => { obj[k] = _nullableValue(cfg[k], false); });

  obj.rnacalibrate = {
    calibration_variant: cfg.rnacalibrate.calibration_variant,
    k: cfg.rnacalibrate.k,
    max_target_length: cfg.rnacalibrate.max_target_length,
    randomize_targets: cfg.rnacalibrate.randomize_targets,
  };

  const distForced = window.isDistributionForced(cfg);
  obj.rnahybrid = { species: cfg.rnahybrid.species };
  window.RNAHYBRID_NULLABLE_KEYS.forEach((k) => {
    obj.rnahybrid[k] = _nullableValue(cfg.rnahybrid[k], k === "distribution" && distForced);
  });

  // Helix values are written even while the form greys them out (model !== "B"), so
  // toggling the model never destroys the user's tuning; the backend warns they are inert.
  const _block = (state, spec) => Object.fromEntries(
    spec.map(([k, kind]) => [
      k,
      kind === "n" ? _nullableValue(state[k], false)
      : kind === "c" ? window.INTARNA_OUTPUT_COLUMNS_DEFAULT
      : state[k],
    ])
  );

  const it = cfg.intarna;
  const B = window.INTARNA_BLOCKS;
  obj.intarna = {
    accessibility_variant: it.accessibility_variant,
    prediction_mode: it.prediction_mode,
    model: it.model,
    energy_set: it.energy_set,
    ..._block(it, B.top),
    helix:         _block(it.helix,         B.helix),
    seed:          _block(it.seed,          B.seed),
    accessibility: _block(it.accessibility, B.accessibility),
    output:        _block(it.output,        B.output),
  };

  // `type` alone enables plots (Snakefile: `plots_enabled = plots_type is not None`);
  // every other key is read only inside rule build_plots. The rest of the block is
  // therefore inert under `type: null`, so always write it — disabling plots must not
  // destroy the user's basesize/thresholds/locus lists.
  const pl = cfg.plots;
  const allTypes = window.PLOT_TYPES.map((p) => p.value);
  const isAll = pl.types.length > 0 && allTypes.every((t) => pl.types.includes(t));
  obj.plots = {
    type: pl.types.length === 0 ? null : (isAll ? "all" : pl.types.join(",")),
    basesize: pl.basesize,
    pvalue_threshold: pl.pvalue_threshold,
    per_mirna_top_n: pl.per_mirna_top_n,
    locus: _pairListToArray(pl.locus),
    gene: _pairListToArray(pl.gene),
    protein: _pairListToArray(pl.protein),
  };

  obj.results_dir = cfg.results_dir;
  return obj;
};

// ---- inverse: parsed config dict -> form state (start from INITIAL_CONFIG defaults) ----

function _toNullableScalar(v, dflt) {
  return v === null || v === undefined ? { set: false, value: dflt } : { set: true, value: v };
}
function _toNullablePair(v, da, db) {
  if (v === null || v === undefined) return { set: false, a: da, b: db };
  const [a, b] = String(v).split(",");
  return { set: true, a: Number(a), b: Number(b) };
}
// PyYAML (YAML 1.1) turns bare on/off into booleans; coerce variant fields back to strings.
function _variantStr(v, dflt) {
  if (v === true) return "on";
  if (v === false) return "off";
  return typeof v === "string" ? v : dflt;
}
function _arrayToPairList(items, prefix) {
  if (!Array.isArray(items)) return [];
  return items.map((s, i) => {
    const [value, mirna] = String(s).split(",");
    return { id: prefix + i + Date.now().toString(36), value: value || "", mirna: mirna || "" };
  });
}

window.hydrateConfig = function hydrateConfig(raw) {
  const cfg = JSON.parse(JSON.stringify(window.INITIAL_CONFIG));
  if (!raw || typeof raw !== "object") return cfg;

  if (raw.queries && typeof raw.queries === "object") {
    cfg.queries = Object.entries(raw.queries).map(([key, path], i) =>
      ({ id: "q" + i + Date.now().toString(36), key, path: path || "" }));
  }
  if (raw.targets && typeof raw.targets === "object") {
    cfg.targets = Object.entries(raw.targets).map(([key, path], i) =>
      ({ id: "t" + i + Date.now().toString(36), key, path: path || "" }));
  }

  if (typeof raw.threads === "number") cfg.threads = raw.threads;
  cfg.max_suboptimal_hits = _toNullableScalar(raw.max_suboptimal_hits, cfg.max_suboptimal_hits.value);
  cfg.max_total_energy = _toNullableScalar(raw.max_total_energy, cfg.max_total_energy.value);
  cfg.seed = _toNullablePair(raw.seed, cfg.seed.a, cfg.seed.b);

  if (raw.rnacalibrate) {
    const rc = raw.rnacalibrate;
    cfg.rnacalibrate.calibration_variant = _variantStr(rc.calibration_variant, cfg.rnacalibrate.calibration_variant);
    if (rc.k != null) cfg.rnacalibrate.k = rc.k;
    if (rc.max_target_length != null) cfg.rnacalibrate.max_target_length = rc.max_target_length;
    if (rc.randomize_targets != null) cfg.rnacalibrate.randomize_targets = rc.randomize_targets;
  }

  if (raw.rnahybrid) {
    const rh = raw.rnahybrid;
    if (rh.species != null) cfg.rnahybrid.species = rh.species;
    cfg.rnahybrid.max_internal_loop = _toNullableScalar(rh.max_internal_loop, cfg.rnahybrid.max_internal_loop.value);
    cfg.rnahybrid.max_bulge_loop = _toNullableScalar(rh.max_bulge_loop, cfg.rnahybrid.max_bulge_loop.value);
    cfg.rnahybrid.pvalue_threshold = _toNullableScalar(rh.pvalue_threshold, cfg.rnahybrid.pvalue_threshold.value);
    cfg.rnahybrid.distribution = _toNullablePair(rh.distribution, cfg.rnahybrid.distribution.a, cfg.rnahybrid.distribution.b);
  }

  if (raw.intarna) {
    const ri = raw.intarna, ci = cfg.intarna;
    ci.accessibility_variant = _variantStr(ri.accessibility_variant, ci.accessibility_variant);
    if (ri.prediction_mode != null) ci.prediction_mode = ri.prediction_mode;
    if (ri.model != null) ci.model = ri.model;
    if (ri.energy_set != null) ci.energy_set = ri.energy_set;
    // "c" (columns) is a constant the form never edits, so it is not hydrated back.
    const _hydrateBlock = (rawBlock, state, spec) => {
      if (!rawBlock) return;
      spec.forEach(([k, kind]) => {
        if (kind === "n") state[k] = _toNullableScalar(rawBlock[k], state[k].value);
        else if (kind === "p" && rawBlock[k] != null) state[k] = rawBlock[k];
      });
    };

    const B = window.INTARNA_BLOCKS;
    _hydrateBlock(ri,               ci,               B.top);
    _hydrateBlock(ri.helix,         ci.helix,         B.helix);
    _hydrateBlock(ri.seed,          ci.seed,          B.seed);
    _hydrateBlock(ri.accessibility, ci.accessibility, B.accessibility);
    _hydrateBlock(ri.output,        ci.output,        B.output);
  }

  if (raw.plots) {
    const rp = raw.plots;
    // Read the surrounding settings whether or not `type` is null — see the writer.
    const all = window.PLOT_TYPES.map((p) => p.value);
    cfg.plots.types = rp.type == null ? []
      : rp.type === "all" ? all
      : String(rp.type).split(",").map((s) => s.trim()).filter(Boolean);
    if (rp.basesize != null) cfg.plots.basesize = rp.basesize;
    if (rp.pvalue_threshold != null) cfg.plots.pvalue_threshold = rp.pvalue_threshold;
    if (rp.per_mirna_top_n != null) cfg.plots.per_mirna_top_n = rp.per_mirna_top_n;
    if (rp.locus != null) cfg.plots.locus = _arrayToPairList(rp.locus, "l");
    if (rp.gene != null) cfg.plots.gene = _arrayToPairList(rp.gene, "g");
    if (rp.protein != null) cfg.plots.protein = _arrayToPairList(rp.protein, "p");
  }

  if (raw.results_dir != null) cfg.results_dir = raw.results_dir;
  return cfg;
};
