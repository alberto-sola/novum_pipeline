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

// One nullable RNAhybrid field → one IR line. Empty inputs and the
// distribution-forced-null case all collapse to `<key>: null`.
function nullableLine(key, field, forced) {
  if (forced || !field.set) {
    return { t: "nested", indent: 2, k: key, v: null, kind: "null" };
  }
  if (key === "seed" || key === "distribution") {
    const { a, b } = field;
    if (a === "" || b === "" || a == null || b == null) {
      return { t: "nested", indent: 2, k: key, v: null, kind: "null" };
    }
    return { t: "nested", indent: 2, k: key, v: `${a},${b}`, kind: "str" };
  }
  if (field.value === "" || field.value === null || field.value === undefined) {
    return { t: "nested", indent: 2, k: key, v: null, kind: "null" };
  }
  return { t: "nested", indent: 2, k: key, v: field.value, kind: "num" };
}

window.buildYAML = function buildYAML(cfg) {
  const lines = [];

  pushMappingSection(lines, "queries", cfg.queries, "unnamed_query");
  lines.push({ t: "blank" });

  pushMappingSection(lines, "targets", cfg.targets, "unnamed_target");
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnacalibrate" });
  lines.push({ t: "nested", indent: 2, k: "mode",              v: cfg.rnacalibrate.mode,              kind: "str"  });
  lines.push({ t: "nested", indent: 2, k: "k",                 v: cfg.rnacalibrate.k,                 kind: "num"  });
  lines.push({ t: "nested", indent: 2, k: "max_target_length", v: cfg.rnacalibrate.max_target_length, kind: "num"  });
  lines.push({ t: "nested", indent: 2, k: "randomize_targets", v: cfg.rnacalibrate.randomize_targets, kind: "bool" });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnahybrid" });
  lines.push({ t: "nested", indent: 2, k: "threads", v: cfg.rnahybrid.threads, kind: "num" });
  lines.push({ t: "nested", indent: 2, k: "species", v: cfg.rnahybrid.species, kind: "str" });
  const distForced = window.isDistributionForced(cfg);
  window.RNAHYBRID_NULLABLE_KEYS.forEach((key) => {
    lines.push(nullableLine(key, cfg.rnahybrid[key], key === "distribution" && distForced));
  });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "build_plots" });
  const bp = cfg.build_plots;
  if (bp.types.length === 0) {
    // Skip-plotting state: emit only `type: null`. Surrounding fields are
    // intentionally omitted so disabled runs produce a minimal config block.
    lines.push({ t: "nested", indent: 2, k: "type", v: null, kind: "null" });
  } else {
    const allTypes = window.PLOT_TYPES.map((p) => p.value);
    const isAll = allTypes.every((t) => bp.types.includes(t));
    const typeVal = isAll ? "all" : bp.types.join(",");
    lines.push({ t: "nested", indent: 2, k: "type",     v: typeVal,     kind: "str" });
    lines.push({ t: "nested", indent: 2, k: "basesize", v: bp.basesize, kind: "num" });
    lines.push({ t: "nested", indent: 2, k: "pvalue",   v: bp.pvalue,   kind: "num" });

    lines.push({ t: "nested", indent: 2, k: "locus", v: "", kind: "bare" });
    const validLoci = bp.locus
      .map((l) => ({ v: (l.value || "").trim(), m: (l.mirna || "").trim() }))
      .filter((l) => l.v);
    if (validLoci.length === 0) {
      lines.push({ t: "empty_list", indent: 4 });
    } else {
      validLoci.forEach((l) => lines.push({ t: "listitem", indent: 4, v: l.m ? `${l.v},${l.m}` : l.v, kind: "str" }));
    }

    lines.push({ t: "nested", indent: 2, k: "gene", v: "", kind: "bare" });
    const validGenes = bp.gene
      .map((g) => ({ v: (g.value || "").trim(), m: (g.mirna || "").trim() }))
      .filter((g) => g.v);
    if (validGenes.length === 0) {
      lines.push({ t: "empty_list", indent: 4 });
    } else {
      validGenes.forEach((g) => lines.push({ t: "listitem", indent: 4, v: g.m ? `${g.v},${g.m}` : g.v, kind: "str" }));
    }

    lines.push({ t: "nested", indent: 2, k: "protein", v: "", kind: "bare" });
    const validProteins = bp.protein
      .map((p) => ({ v: (p.value || "").trim(), m: (p.mirna || "").trim() }))
      .filter((p) => p.v);
    if (validProteins.length === 0) {
      lines.push({ t: "empty_list", indent: 4 });
    } else {
      validProteins.forEach((p) => lines.push({ t: "listitem", indent: 4, v: p.m ? `${p.v},${p.m}` : p.v, kind: "str" }));
    }
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
    if (L.t === "raw") return `<span>${esc(L.s)}</span>`;
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
    if (L.t === "raw") return L.s;
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
