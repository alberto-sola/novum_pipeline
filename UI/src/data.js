// Novum Pipeline — reference data
window.SPECIES_OPTIONS = [
  { value: "3utr_human", label: "Human (3′UTR)" },
  { value: "3utr_fly",   label: "Fly (3′UTR)" },
  { value: "3utr_worm",  label: "Worm (3′UTR)" }
];

window.CALIBRATE_MODES = [
  { value: "calibrated",   label: "Calibrated" },
  { value: "uncalibrated", label: "Uncalibrated" },
  { value: "both",         label: "Both" }
];

window.INTARNA_ACC_MODES = [
  { value: "on",   label: "On"   },
  { value: "off",  label: "Off"  },
  { value: "both", label: "Both" }
];

window.INTARNA_PREDICTION_MODES = [
  { value: "H", label: "H · heuristic" },
  { value: "M", label: "M · exact"     },
  { value: "S", label: "S · seed-only" }
];

window.INTARNA_MODELS = [
  { value: "X", label: "X · seed-extension MFE" },
  { value: "S", label: "S · legacy MFE"          },
  { value: "B", label: "B · helix-block"         },
  { value: "P", label: "P · ensemble"            }
];

window.INTARNA_OVERLAPS = [
  { value: "N", label: "N · none"          },
  { value: "T", label: "T · target only"   },
  { value: "Q", label: "Q · query only"    },
  { value: "B", label: "B · both"          }
];

window.PLOT_TYPES = [
  { value: "pvalue_distribution", label: "p-value distribution" },
  { value: "position_pvalue",     label: "position × p-value" },
  { value: "position_energy",     label: "position × energy" },
  { value: "volcano",             label: "volcano" },
  { value: "pvalue_ecdf",         label: "p-value ECDF" },
  { value: "calibration_delta",   label: "calibration Δ" },
  { value: "per_mirna",           label: "per-miRNA" },
  { value: "position_density",    label: "position density" },
  { value: "seed_class",          label: "seed class" }
];

// Nullable-field key lists per card. Used by section renderers (and the Hero
// "set params" counter) so a single source of truth drives what's rendered
// vs. what's just present in INITIAL_CONFIG.
window.SHARED_NULLABLE_KEYS    = ["max_suboptimal_hits", "seed", "max_total_energy"];
window.RNAHYBRID_NULLABLE_KEYS = ["max_internal_loop", "max_bulge_loop", "pvalue_threshold", "distribution"];
window.INTARNA_SEED_NULLABLE_KEYS = [
  "length", "query_range", "target_range", "max_unpaired_bases",
  "max_energy", "max_hybrid_energy", "min_unpaired_probability"
];
window.INTARNA_ACC_NULLABLE_KEYS    = ["window", "max_bp_span"];
window.INTARNA_OUTPUT_NULLABLE_KEYS = ["max_delta_energy", "min_unpaired_probability"];

// `distribution` is supplied by RNAcalibrate whenever calibration is part of the run,
// so the RNAhybrid form forces it null in any variant other than uncalibrated.
window.isDistributionForced = (cfg) => cfg.rnacalibrate.calibration_variant !== "uncalibrated";

// The three IntaRNA seed sub-fields below derive from the top-level `seed` at runtime
// (intarna.py:_add_seed_flags), so the UI hard-locks them whenever the shared seed is set.
window.isIntarnaSeedDerived = (cfg) => cfg.seed.set;

// Default IntaRNA --outCsvCols string from Config/config.yaml. The form does not
// expose `intarna.output.columns`; the YAML writer emits this constant.
window.INTARNA_OUTPUT_COLUMNS_DEFAULT =
  "id1,id2,start1,end1,start2,end2,subseqDP,hybridDP,E,E_hybrid,ED1,ED2,Pu1,Pu2,seedStart1,seedEnd1,seedE,seedStart2,seedEnd2";

// canonical initial state — queries/targets populated from the launcher's dropdowns.
// Plots are toggled by the size of `plots.types`: an empty list serializes to
// `type: null` (skip), any selection emits the full plot block.
window.INITIAL_CONFIG = {
  queries: [],
  targets: [],
  threads: 16,
  max_suboptimal_hits: { set: false, value: 1 },
  seed:                { set: false, a: 2, b: 7 },
  max_total_energy:    { set: false, value: -18 },
  rnacalibrate: {
    calibration_variant: "calibrated",
    k: 10000,
    max_target_length: 50000,
    randomize_targets: true
  },
  rnahybrid: {
    species: "3utr_human",
    max_internal_loop: { set: false, value: 1 },
    max_bulge_loop:    { set: false, value: 1 },
    pvalue_threshold:  { set: false, value: 0.05 },
    distribution:      { set: false, a: 0, b: 1 }
  },
  intarna: {
    accessibility_variant: "on",
    prediction_mode: "H",
    model: "S",
    max_interaction_length: 0,
    max_loop_size: 10,
    seed: {
      enabled: false,
      length:                   { set: false, value: 7 },
      max_energy:               { set: false, value: 0 },
      max_hybrid_energy:        { set: false, value: 999 },
      min_unpaired_probability: { set: false, value: 0 },
      forbid_gu: false,
      forbid_gu_at_ends: false,
      query_range:              { set: false, value: "" },
      target_range:             { set: false, value: "" },
      max_unpaired_bases:       { set: false, value: 0 },
      report_best_only: false
    },
    accessibility: {
      window:      { set: false, value: 150 },
      max_bp_span: { set: false, value: 100 },
      forbid_lonely_pairs: false,
      forbid_gu_at_ends: false
    },
    output: {
      max_delta_energy:         { set: false, value: 100 },
      overlap: "B",
      min_unpaired_probability: { set: false, value: 0 },
      forbid_lonely_pairs: false,
      forbid_gu_at_ends: false
    }
  },
  plots: {
    types: window.PLOT_TYPES.map((p) => p.value),
    basesize: 12,
    pvalue_threshold: 0.01,
    per_mirna_top_n: 12,
    locus: [],
    gene: [],
    protein: []
  },
  results_dir: "Data/Results/"
};
