// Novum Pipeline — YAML serializer (matches Config/config.yaml schema)
window.buildYAML = function buildYAML(cfg) {
  const lines = [];

  lines.push({ t: "k", k: "queries" });
  if (!cfg.queries.length) {
    lines.push({ t: "raw", s: "  {}" });
  } else {
    cfg.queries.forEach((q) => {
      const k = (q.key || "").trim() || "unnamed_query";
      lines.push({ t: "nested", indent: 2, k, v: q.path || "", kind: "str" });
    });
  }
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "targets" });
  if (!cfg.targets.length) {
    lines.push({ t: "raw", s: "  {}" });
  } else {
    cfg.targets.forEach((t) => {
      const k = (t.key || "").trim() || "unnamed_target";
      lines.push({ t: "nested", indent: 2, k, v: t.path || "", kind: "str" });
    });
  }
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnacalibrate" });
  lines.push({ t: "nested", indent: 2, k: "mode", v: cfg.rnacalibrate.mode, kind: "str" });
  lines.push({ t: "nested", indent: 2, k: "k", v: cfg.rnacalibrate.k, kind: "num" });
  lines.push({ t: "nested", indent: 2, k: "max_target_length", v: cfg.rnacalibrate.max_target_length, kind: "num" });
  lines.push({ t: "nested", indent: 2, k: "randomize_targets", v: cfg.rnacalibrate.randomize_targets, kind: "bool" });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "rnahybrid" });
  lines.push({ t: "nested", indent: 2, k: "threads", v: cfg.rnahybrid.threads, kind: "num" });
  lines.push({ t: "nested", indent: 2, k: "species", v: cfg.rnahybrid.species, kind: "str" });
  ["hits", "u", "v", "energy", "pvalue", "seed", "distribution"].forEach((key) => {
    const f = cfg.rnahybrid[key];
    const forcedNull = (key === "distribution" && cfg.rnacalibrate.mode !== "uncalibrated");
    if (forcedNull || !f.set) {
      lines.push({ t: "nested", indent: 2, k: key, v: null, kind: "null" });
      return;
    }
    if (key === "seed" || key === "distribution") {
      const a = f.a, b = f.b;
      if (a === "" || b === "" || a == null || b == null) {
        lines.push({ t: "nested", indent: 2, k: key, v: null, kind: "null" });
      } else {
        lines.push({ t: "nested", indent: 2, k: key, v: `${a},${b}`, kind: "str" });
      }
    } else {
      if (f.value === "" || f.value === null || f.value === undefined) {
        lines.push({ t: "nested", indent: 2, k: key, v: null, kind: "null" });
      } else {
        lines.push({ t: "nested", indent: 2, k: key, v: f.value, kind: "num" });
      }
    }
  });
  lines.push({ t: "blank" });

  lines.push({ t: "raw", s: "tidy_rnahybrid: {}" });
  lines.push({ t: "blank" });

  lines.push({ t: "raw", s: "annotate_rnahybrid: {}" });
  lines.push({ t: "blank" });

  lines.push({ t: "k", k: "build_plots" });
  const bp = cfg.build_plots;
  if (!bp.enabled) {
    lines.push({ t: "nested", indent: 2, k: "type", v: null, kind: "null" });
  } else {
    const allTypes = window.PLOT_TYPES.map((p) => p.value);
    const isAll = allTypes.every((t) => bp.types.includes(t));
    const typeVal = bp.types.length === 0 ? null : (isAll ? "all" : bp.types.join(","));
    lines.push({ t: "nested", indent: 2, k: "type", v: typeVal, kind: typeVal == null ? "null" : "str" });
    lines.push({ t: "nested", indent: 2, k: "basesize", v: bp.basesize, kind: "num" });
    lines.push({ t: "nested", indent: 2, k: "pvalue", v: bp.pvalue, kind: "num" });

    lines.push({ t: "nested", indent: 2, k: "locus", v: "", kind: "bare" });
    const validLoci = bp.locus.map((l) => (l.value || "").trim()).filter(Boolean);
    if (validLoci.length === 0) {
      lines.push({ t: "raw", s: "    []" });
    } else {
      validLoci.forEach((v) => {
        lines.push({ t: "listitem", indent: 4, v, kind: "str" });
      });
    }

    lines.push({ t: "nested", indent: 2, k: "gene", v: "", kind: "bare" });
    const validGenes = bp.gene
      .map((g) => ({ v: (g.value || "").trim(), m: (g.mirna || "").trim() }))
      .filter((g) => g.v);
    if (validGenes.length === 0) {
      lines.push({ t: "raw", s: "    []" });
    } else {
      validGenes.forEach((g) => {
        lines.push({ t: "listitem", indent: 4, v: g.m ? `${g.v},${g.m}` : g.v, kind: "str" });
      });
    }

    // protein UI deferred — emit an empty key to mirror Config/config.yaml.
    lines.push({ t: "nested", indent: 2, k: "protein", v: "", kind: "bare" });
  }
  lines.push({ t: "blank" });

  lines.push({ t: "kv", k: "results_dir", v: cfg.results_dir, kind: "str" });

  return lines;
};

const esc = (s) => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

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
    if (L.t === "nested") {
      const head = `${" ".repeat(L.indent)}<span class="y-key">${esc(L.k)}</span>:`;
      return L.kind === "bare" ? head : `${head} ${val(L.kind, L.v)}`;
    }
    if (L.t === "listitem") {
      return `${" ".repeat(L.indent)}- ${val(L.kind, L.v)}`;
    }
    return "";
  }).join("\n");
};

window.renderYAMLPlain = function renderYAMLPlain(lines) {
  const v = (kind, v) => {
    if (kind === "bare") return "";
    if (kind === "null" || v === null || v === undefined) return "null";
    if (kind === "bool") return v ? "true" : "false";
    return String(v);
  };
  return lines.map((L) => {
    if (L.t === "blank") return "";
    if (L.t === "raw") return L.s;
    if (L.t === "k") return `${L.k}:`;
    if (L.t === "kv") return `${L.k}: ${v(L.kind, L.v)}`;
    if (L.t === "nested") {
      const head = `${" ".repeat(L.indent)}${L.k}:`;
      return L.kind === "bare" ? head : `${head} ${v(L.kind, L.v)}`;
    }
    if (L.t === "listitem") {
      return `${" ".repeat(L.indent)}- ${v(L.kind, L.v)}`;
    }
    return "";
  }).join("\n");
};
