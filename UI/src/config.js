// Novum Pipeline — config <-> form-state bridge (mirrors yaml.js buildYAML schema)

// form-state nullable field -> its plain value (or null). Delegates to yaml.js's
// resolveNullable so the null-collapsing rules live in exactly one place.
function _nullableValue(field, forced) {
  return window.resolveNullable(field, forced).v;
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

  const it = cfg.intarna;
  obj.intarna = {
    accessibility_variant: it.accessibility_variant,
    prediction_mode: it.prediction_mode,
    model: it.model,
    energy_set: it.energy_set,
    max_interaction_length: it.max_interaction_length,
    max_loop_size: it.max_loop_size,
    helix: {
      // Written even while greyed out (model !== "B"), so toggling the model to compare
      // arms never destroys the user's tuning; intarna.py warns that they are inert.
      ...Object.fromEntries(
        window.INTARNA_HELIX_NULLABLE_KEYS.map((k) => [k, _nullableValue(it.helix[k], false)])
      ),
      full_energy: it.helix.full_energy,
    },
    seed: {
      max_energy: _nullableValue(it.seed.max_energy, false),
      max_hybrid_energy: _nullableValue(it.seed.max_hybrid_energy, false),
      min_unpaired_probability: _nullableValue(it.seed.min_unpaired_probability, false),
      forbid_gu: it.seed.forbid_gu,
      forbid_gu_at_ends: it.seed.forbid_gu_at_ends,
      target_range: _nullableValue(it.seed.target_range, false),
      max_unpaired_bases: _nullableValue(it.seed.max_unpaired_bases, false),
      report_best_only: it.seed.report_best_only,
    },
    accessibility: {
      window: _nullableValue(it.accessibility.window, false),
      max_bp_span: _nullableValue(it.accessibility.max_bp_span, false),
      forbid_lonely_pairs: it.accessibility.forbid_lonely_pairs,
      forbid_gu_at_ends: it.accessibility.forbid_gu_at_ends,
    },
    output: {
      max_delta_energy: _nullableValue(it.output.max_delta_energy, false),
      overlap: it.output.overlap,
      min_unpaired_probability: _nullableValue(it.output.min_unpaired_probability, false),
      forbid_lonely_pairs: it.output.forbid_lonely_pairs,
      forbid_gu_at_ends: it.output.forbid_gu_at_ends,
      columns: window.INTARNA_OUTPUT_COLUMNS_DEFAULT,
    },
  };

  const pl = cfg.plots;
  if (pl.types.length === 0) {
    obj.plots = { type: null };
  } else {
    const allTypes = window.PLOT_TYPES.map((p) => p.value);
    const isAll = allTypes.every((t) => pl.types.includes(t));
    obj.plots = {
      type: isAll ? "all" : pl.types.join(","),
      basesize: pl.basesize,
      pvalue_threshold: pl.pvalue_threshold,
      per_mirna_top_n: pl.per_mirna_top_n,
      locus: _pairListToArray(pl.locus),
      gene: _pairListToArray(pl.gene),
      protein: _pairListToArray(pl.protein),
    };
  }

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
    if (ri.max_interaction_length != null) ci.max_interaction_length = ri.max_interaction_length;
    if (ri.max_loop_size != null) ci.max_loop_size = ri.max_loop_size;
    if (ri.energy_set != null) ci.energy_set = ri.energy_set;
    if (ri.helix) {
      window.INTARNA_HELIX_NULLABLE_KEYS.forEach((k) => {
        ci.helix[k] = _toNullableScalar(ri.helix[k], ci.helix[k].value);
      });
      if (ri.helix.full_energy != null) ci.helix.full_energy = ri.helix.full_energy;
    }
    if (ri.seed) {
      ci.seed.max_energy = _toNullableScalar(ri.seed.max_energy, ci.seed.max_energy.value);
      ci.seed.max_hybrid_energy = _toNullableScalar(ri.seed.max_hybrid_energy, ci.seed.max_hybrid_energy.value);
      ci.seed.min_unpaired_probability = _toNullableScalar(ri.seed.min_unpaired_probability, ci.seed.min_unpaired_probability.value);
      if (ri.seed.forbid_gu != null) ci.seed.forbid_gu = ri.seed.forbid_gu;
      if (ri.seed.forbid_gu_at_ends != null) ci.seed.forbid_gu_at_ends = ri.seed.forbid_gu_at_ends;
      ci.seed.target_range = _toNullableScalar(ri.seed.target_range, ci.seed.target_range.value);
      ci.seed.max_unpaired_bases = _toNullableScalar(ri.seed.max_unpaired_bases, ci.seed.max_unpaired_bases.value);
      if (ri.seed.report_best_only != null) ci.seed.report_best_only = ri.seed.report_best_only;
    }
    if (ri.accessibility) {
      ci.accessibility.window = _toNullableScalar(ri.accessibility.window, ci.accessibility.window.value);
      ci.accessibility.max_bp_span = _toNullableScalar(ri.accessibility.max_bp_span, ci.accessibility.max_bp_span.value);
      if (ri.accessibility.forbid_lonely_pairs != null) ci.accessibility.forbid_lonely_pairs = ri.accessibility.forbid_lonely_pairs;
      if (ri.accessibility.forbid_gu_at_ends != null) ci.accessibility.forbid_gu_at_ends = ri.accessibility.forbid_gu_at_ends;
    }
    if (ri.output) {
      ci.output.max_delta_energy = _toNullableScalar(ri.output.max_delta_energy, ci.output.max_delta_energy.value);
      if (ri.output.overlap != null) ci.output.overlap = ri.output.overlap;
      ci.output.min_unpaired_probability = _toNullableScalar(ri.output.min_unpaired_probability, ci.output.min_unpaired_probability.value);
      if (ri.output.forbid_lonely_pairs != null) ci.output.forbid_lonely_pairs = ri.output.forbid_lonely_pairs;
      if (ri.output.forbid_gu_at_ends != null) ci.output.forbid_gu_at_ends = ri.output.forbid_gu_at_ends;
    }
  }

  if (raw.plots) {
    const rp = raw.plots;
    if (rp.type == null) {
      cfg.plots.types = [];
    } else {
      const all = window.PLOT_TYPES.map((p) => p.value);
      cfg.plots.types = rp.type === "all" ? all : String(rp.type).split(",").map((s) => s.trim()).filter(Boolean);
      if (rp.basesize != null) cfg.plots.basesize = rp.basesize;
      if (rp.pvalue_threshold != null) cfg.plots.pvalue_threshold = rp.pvalue_threshold;
      if (rp.per_mirna_top_n != null) cfg.plots.per_mirna_top_n = rp.per_mirna_top_n;
      cfg.plots.locus = _arrayToPairList(rp.locus, "l");
      cfg.plots.gene = _arrayToPairList(rp.gene, "g");
      cfg.plots.protein = _arrayToPairList(rp.protein, "p");
    }
  }

  if (raw.results_dir != null) cfg.results_dir = raw.results_dir;
  return cfg;
};
