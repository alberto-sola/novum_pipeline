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

// canonical initial state — queries/targets populated from the launcher's dropdowns
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
    enabled: true,
    types: ["pvalue_distribution", "position_pvalue", "position_energy"], // "all" when all selected
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
