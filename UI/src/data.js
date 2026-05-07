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

window.RNAHYBRID_NULLABLE_KEYS = ["hits", "u", "v", "energy", "pvalue", "seed", "distribution"];

// `distribution` is supplied by RNAcalibrate whenever calibration is part of the run,
// so the RNAhybrid form forces it null in any mode other than uncalibrated.
window.isDistributionForced = (cfg) => cfg.rnacalibrate.mode !== "uncalibrated";

// canonical initial state — queries/targets populated from the launcher's dropdowns.
// Plots are toggled by the size of `build_plots.types`: an empty list serializes to
// `type: null` (skip), any selection emits the full plot block.
window.INITIAL_CONFIG = {
  queries: [],
  targets: [],
  rnacalibrate: {
    mode: "calibrated",
    k: 10000,
    max_target_length: 50000,
    randomize_targets: true
  },
  rnahybrid: {
    threads: 8,
    species: "3utr_human",
    hits:         { set: false, value: 10 },
    u:            { set: false, value: 1 },
    v:            { set: false, value: 1 },
    energy:       { set: false, value: -18 },
    pvalue:       { set: false, value: 0.05 },
    seed:         { set: false, a: 2, b: 7 },
    distribution: { set: false, a: 0, b: 1 }
  },
  build_plots: {
    types: window.PLOT_TYPES.map((p) => p.value),
    basesize: 12,
    pvalue: 0.01,
    per_mirna_top_n: 12,
    locus: [],
    gene: [],
    protein: []
  },
  results_dir: "Data/Results/"
};
