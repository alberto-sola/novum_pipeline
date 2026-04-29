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
  { value: "position_energy",     label: "position × energy" }
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
    mode: "both",
    k: 10000,
    max_target_length: 50000,
    randomize_targets: true
  },
  rnahybrid: {
    threads: 16,
    species: "3utr_human",
    hits:         { set: false, value: 10 },
    u:            { set: false, value: 1 },
    v:            { set: false, value: 1 },
    energy:       { set: true,  value: -18 },
    pvalue:       { set: false, value: 0.05 },
    seed:         { set: false, a: 2, b: 7 },
    distribution: { set: false, a: 0, b: 1 }
  },
  build_plots: {
    types: window.PLOT_TYPES.map((p) => p.value),
    basesize: 12,
    pvalue: 0.01,
    locus: [
      { id: "l1", value: "FE838_RS16060" },
      { id: "l2", value: "FE838_RS16090" },
      { id: "l3", value: "FE838_RS16070" },
      { id: "l4", value: "FE838_RS16075" },
      { id: "l5", value: "FE838_RS16085" },
      { id: "l6", value: "FE838_RS16065" },
      { id: "l7", value: "FE838_RS16080" }
    ],
    gene: [
      { id: "g1", value: "yegH", mirna: "hsa-miR-1226-5p" },
      { id: "g2", value: "rnpA", mirna: "hsa-miR-4747-3p" }
    ]
  },
  results_dir: "Data/Results/"
};
