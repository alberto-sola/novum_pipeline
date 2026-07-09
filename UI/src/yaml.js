// Novum Pipeline — YAML serializer (matches Config/config.yaml schema)

function pushMappingSection(lines, key, items, fallbackKey) {
  lines.push({ t: "k", k: key });
  if (!items.length) {
    lines.push({ t: "empty_map", indent: 2 });
    return;
  }
  items.forEach((item) => {
    const k = (item.key || "").trim() || fallbackKey;
    lines.push({ t: "nested", indent: 2, k, v: item.path || "", kind: "str" });
  });
}

// Resolve a nullable form-state field to its serialized {v, kind}. Pair-shaped
// fields (with `a`/`b`) emit "a,b" strings; scalar fields use typeof to pick str
// vs num; unset/forced/empty fields collapse to a null value.
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
  const kind = typeof field.value === "string" ? "str" : "num";
  return { v: field.value, kind };
}

// One nullable field → one IR line. Generalised over indent so it works for
// top-level shared params, rnahybrid params, and the IntaRNA sub-blocks.
function nullableLineAt(indent, key, field, forced) {
  const { v, kind } = resolveNullable(field, forced);
  return { t: "nested", indent, k: key, v, kind };
}

// One plots pair-list (value + optional miRNA) → a `key:` header then either an
// empty list or one `value` / `value,mirna` item per non-blank entry.
function pushPairList(lines, key, items) {
  lines.push({ t: "nested", indent: 2, k: key, v: "", kind: "bare" });
  const valid = items
    .map((it) => ({ v: (it.value || "").trim(), m: (it.mirna || "").trim() }))
    .filter((it) => it.v);
  if (valid.length === 0) {
    lines.push({ t: "empty_list", indent: 4 });
  } else {
    valid.forEach((it) => lines.push({ t: "listitem", indent: 4, v: it.m ? `${it.v},${it.m}` : it.v, kind: "str" }));
  }
}

window.buildYAML = function buildYAML(cfg) {
  const lines = [];

  pushMappingSection(lines, "queries", cfg.queries, "unnamed_query");
  lines.push({ t: "blank" });

  pushMappingSection(lines, "targets", cfg.targets, "unnamed_target");
  lines.push({ t: "blank" });

  // shared-by-both-arms block
  lines.push({ t: "kv", k: "threads", v: cfg.threads, kind: "num" });
  window.SHARED_NULLABLE_KEYS.forEach((key) => {
    const { v, kind } = resolveNullable(cfg[key], false);
    lines.push({ t: "kv", k: key, v, kind });
  });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnacalibrate" });
  lines.push({ t: "nested", indent: 2, k: "calibration_variant", v: cfg.rnacalibrate.calibration_variant, kind: "str"  });
  lines.push({ t: "nested", indent: 2, k: "k",                   v: cfg.rnacalibrate.k,                   kind: "num"  });
  lines.push({ t: "nested", indent: 2, k: "max_target_length",   v: cfg.rnacalibrate.max_target_length,   kind: "num"  });
  lines.push({ t: "nested", indent: 2, k: "randomize_targets",   v: cfg.rnacalibrate.randomize_targets,   kind: "bool" });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnahybrid" });
  lines.push({ t: "nested", indent: 2, k: "species", v: cfg.rnahybrid.species, kind: "str" });
  const distForced = window.isDistributionForced(cfg);
  window.RNAHYBRID_NULLABLE_KEYS.forEach((key) => {
    lines.push(nullableLineAt(2, key, cfg.rnahybrid[key], key === "distribution" && distForced));
  });
  lines.push({ t: "blank" });

  // intarna block
  const it = cfg.intarna;
  lines.push({ t: "k", k: "intarna" });
  lines.push({ t: "nested", indent: 2, k: "accessibility_variant",  v: it.accessibility_variant,  kind: "str" });
  lines.push({ t: "nested", indent: 2, k: "prediction_mode",        v: it.prediction_mode,        kind: "str" });
  lines.push({ t: "nested", indent: 2, k: "model",                  v: it.model,                  kind: "str" });
  lines.push({ t: "nested", indent: 2, k: "max_interaction_length", v: it.max_interaction_length, kind: "num" });
  lines.push({ t: "nested", indent: 2, k: "max_loop_size",          v: it.max_loop_size,          kind: "num" });
  lines.push({ t: "blank" });

  // seed sub-block. enabled/length/query_range are inherited from the top-level
  // `seed`, so they are not emitted here; the rest are independent overrides.
  lines.push({ t: "nested", indent: 2, k: "seed", v: "", kind: "bare" });
  lines.push(nullableLineAt(4, "max_energy",               it.seed.max_energy,               false));
  lines.push(nullableLineAt(4, "max_hybrid_energy",        it.seed.max_hybrid_energy,        false));
  lines.push(nullableLineAt(4, "min_unpaired_probability", it.seed.min_unpaired_probability, false));
  lines.push({ t: "nested", indent: 4, k: "forbid_gu",         v: it.seed.forbid_gu,         kind: "bool" });
  lines.push({ t: "nested", indent: 4, k: "forbid_gu_at_ends", v: it.seed.forbid_gu_at_ends, kind: "bool" });
  lines.push(nullableLineAt(4, "target_range",       it.seed.target_range,       false));
  lines.push(nullableLineAt(4, "max_unpaired_bases", it.seed.max_unpaired_bases, false));
  lines.push({ t: "nested", indent: 4, k: "report_best_only",  v: it.seed.report_best_only,  kind: "bool" });
  lines.push({ t: "blank" });

  // accessibility sub-block
  lines.push({ t: "nested", indent: 2, k: "accessibility", v: "", kind: "bare" });
  lines.push(nullableLineAt(4, "window",      it.accessibility.window,      false));
  lines.push(nullableLineAt(4, "max_bp_span", it.accessibility.max_bp_span, false));
  lines.push({ t: "nested", indent: 4, k: "forbid_lonely_pairs", v: it.accessibility.forbid_lonely_pairs, kind: "bool" });
  lines.push({ t: "nested", indent: 4, k: "forbid_gu_at_ends",   v: it.accessibility.forbid_gu_at_ends,   kind: "bool" });
  lines.push({ t: "blank" });

  // output sub-block (columns is a constant emitted from INTARNA_OUTPUT_COLUMNS_DEFAULT)
  lines.push({ t: "nested", indent: 2, k: "output", v: "", kind: "bare" });
  lines.push(nullableLineAt(4, "max_delta_energy",         it.output.max_delta_energy,         false));
  lines.push({ t: "nested", indent: 4, k: "overlap", v: it.output.overlap, kind: "str" });
  lines.push(nullableLineAt(4, "min_unpaired_probability", it.output.min_unpaired_probability, false));
  lines.push({ t: "nested", indent: 4, k: "forbid_lonely_pairs", v: it.output.forbid_lonely_pairs, kind: "bool" });
  lines.push({ t: "nested", indent: 4, k: "forbid_gu_at_ends",   v: it.output.forbid_gu_at_ends,   kind: "bool" });
  lines.push({ t: "nested", indent: 4, k: "columns", v: window.INTARNA_OUTPUT_COLUMNS_DEFAULT, kind: "str" });
  lines.push({ t: "blank" });

  // extra_args is a UI escape hatch; always emitted as an empty list.
  lines.push({ t: "nested", indent: 2, k: "extra_args", v: "", kind: "bare" });
  lines.push({ t: "empty_list", indent: 4 });
  lines.push({ t: "blank" });

  // plots block (was build_plots; pvalue → pvalue_threshold)
  lines.push({ t: "k", k: "plots" });
  const pl = cfg.plots;
  if (pl.types.length === 0) {
    // Skip-plotting state: emit only `type: null`. Surrounding fields are
    // intentionally omitted so disabled runs produce a minimal config block.
    lines.push({ t: "nested", indent: 2, k: "type", v: null, kind: "null" });
  } else {
    const allTypes = window.PLOT_TYPES.map((p) => p.value);
    const isAll = allTypes.every((t) => pl.types.includes(t));
    const typeVal = isAll ? "all" : pl.types.join(",");
    lines.push({ t: "nested", indent: 2, k: "type",             v: typeVal,            kind: "str" });
    lines.push({ t: "nested", indent: 2, k: "basesize",         v: pl.basesize,        kind: "num" });
    lines.push({ t: "nested", indent: 2, k: "pvalue_threshold", v: pl.pvalue_threshold, kind: "num" });
    lines.push({ t: "nested", indent: 2, k: "per_mirna_top_n",  v: pl.per_mirna_top_n, kind: "num" });

    pushPairList(lines, "locus", pl.locus);
    pushPairList(lines, "gene", pl.gene);
    pushPairList(lines, "protein", pl.protein);
  }
  lines.push({ t: "blank" });

  lines.push({ t: "kv", k: "results_dir", v: cfg.results_dir, kind: "str" });

  return lines;
};

const esc = (s) => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

const indent = (n) => " ".repeat(n);

window.renderYAML = function renderYAML(lines) {
  const val = (kind, v) => {
    if (kind === "bare") return "";
    if (kind === "null" || v === null || v === undefined) return `<span class="y-null">null</span>`;
    if (kind === "bool") return `<span class="y-bool">${v ? "true" : "false"}</span>`;
    if (kind === "num")  return `<span class="y-num">${esc(v)}</span>`;
    return `<span class="y-str">${esc(v)}</span>`;
  };
  return lines.map((L) => {
    if (L.t === "blank") return "";
    if (L.t === "k") return `<span class="y-key">${esc(L.k)}</span>:`;
    if (L.t === "kv") return `<span class="y-key">${esc(L.k)}</span>: ${val(L.kind, L.v)}`;
    if (L.t === "empty_map")  return `${indent(L.indent)}{}`;
    if (L.t === "empty_list") return `${indent(L.indent)}[]`;
    if (L.t === "nested") {
      const head = `${indent(L.indent)}<span class="y-key">${esc(L.k)}</span>:`;
      return L.kind === "bare" ? head : `${head} ${val(L.kind, L.v)}`;
    }
    if (L.t === "listitem") {
      return `${indent(L.indent)}- ${val(L.kind, L.v)}`;
    }
    return "";
  }).join("\n");
};

window.renderYAMLPlain = function renderYAMLPlain(lines) {
  const v = (kind, val) => {
    if (kind === "bare") return "";
    if (kind === "null" || val === null || val === undefined) return "null";
    if (kind === "bool") return val ? "true" : "false";
    return String(val);
  };
  return lines.map((L) => {
    if (L.t === "blank") return "";
    if (L.t === "k") return `${L.k}:`;
    if (L.t === "kv") return `${L.k}: ${v(L.kind, L.v)}`;
    if (L.t === "empty_map")  return `${indent(L.indent)}{}`;
    if (L.t === "empty_list") return `${indent(L.indent)}[]`;
    if (L.t === "nested") {
      const head = `${indent(L.indent)}${L.k}:`;
      return L.kind === "bare" ? head : `${head} ${v(L.kind, L.v)}`;
    }
    if (L.t === "listitem") {
      return `${indent(L.indent)}- ${v(L.kind, L.v)}`;
    }
    return "";
  }).join("\n");
};
