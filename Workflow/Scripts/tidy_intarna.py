import pandas as pd

from _common import iter_fasta_records, parse_header_id, ensure_parent
from _length import add_length_corrected, corrected_name


OUTPUT_COLUMNS = [
    "miRNA", "Gene", "E", "E_hybrid", corrected_name("E_hybrid"), "ED1", "ED2", "Pu1", "Pu2",
    "Gene_length", "Start1", "End1", "Start2", "End2", "Position",
    "subseqDP", "hybridDP",
    "seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2",
]

# Columns that must be present in IntaRNA CSV output for this script to work.
# E_hybrid is the gate and sort key, so a trimmed intarna.output.columns must fail here
# rather than silently produce an ungated file.
REQUIRED_COLS = {"id1", "id2", "start1", "end1", "start2", "end2", "E", "E_hybrid"}

# Rows per read batch: ~100 MB, already below the surviving set's own footprint (the one term
# chunking cannot remove), and throughput is flat from 100k to 1M rows — a pure memory knob.
READ_CHUNK_ROWS = 250_000

# One pair. Spelled once because the cap is applied twice — once per batch, once globally —
# and the two must group identically or the pre-cap stops being a no-op on the result.
CAP_KEYS = ["miRNA", "Gene"]

# IntaRNA writes uppercase NAN in the seed columns under --noSeed (top-level seed: null),
# and that spelling is outside pandas' default NA set — the five columns would land as
# Python strings and cost 34% of the frame. Declaring it types them as float, so they also
# reach the enhance report as "NA" rather than the literal string "NAN".
INTARNA_NA_VALUES = ["NAN"]


#----- Streams the target FASTA once and returns a {Gene → sequence length} map for downstream position normalization -----#
def _parse_gene_lengths(fasta_path):
    return {parse_header_id(header): len(seq) for header, seq in iter_fasta_records(fasta_path)}


#----- Column guard, run once against the header rather than per batch. Pu1/Pu2 stay out
#      of REQUIRED_COLS: trimming intarna.output.columns is legitimate. But a floor
#      against a column that is not there would silently admit every row -----#
def _require_intarna_columns(columns, floors):
    missing = REQUIRED_COLS - set(columns)
    if missing:
        raise ValueError(f"IntaRNA CSV is missing required columns: {sorted(missing)}")

    for column, floor in floors:
        if floor is not None and column not in columns:
            raise ValueError(
                f"IntaRNA CSV is missing {column}, required by the accessibility floor "
                f"({floor}). Keep {column} in intarna.output.columns or clear the floor."
            )


#----- Everything row-independent — rename, validate, gate, floor, annotate — so it can run
#      per batch and let all but the survivors go. The sort and the per-pair cap need the
#      whole surviving set and stay in tidy_intarna below -----#
def _prepare_chunk(chunk, gene_lengths, max_hybrid_energy, floors):
    chunk = chunk.rename(columns={"id1": "Gene", "id2": "miRNA",
                                  "start1": "Start1", "end1": "End1",
                                  "start2": "Start2", "end2": "End2"})

    # Validated on the DISTINCT genes, before any filtering, so a target FASTA that does
    # not cover the output still fails even when the offending rows would have been gated
    # away. Row counts are only computed on the error path.
    genes = chunk["Gene"].unique()
    unknown = [gene for gene in genes if gene not in gene_lengths]
    if unknown:
        count = int(chunk["Gene"].isin(unknown).sum())
        raise ValueError(f"Gene_length missing for {count} record(s): {unknown[:5]}")

    nonpositive = [gene for gene in genes if gene_lengths[gene] <= 0]
    if nonpositive:
        raise ValueError(f"Non-positive Gene_length in target FASTA for: {nonpositive}")

    # Gate on hybridization energy — the quantity RNAhybrid's -e filters and the scale the
    # literature's -18 threshold is stated on. --outMaxE cannot express this under acc=C,
    # where it bounds E_hybrid+ED1+ED2 instead, so this is the authoritative gate on both
    # variants. On wo_accessibility it is a no-op: the tool already applied the same bound.
    if max_hybrid_energy is not None:
        chunk = chunk[chunk["E_hybrid"] <= float(max_hybrid_energy)]

    # Before the cap, like the energy gate: capping first collapses each pair to its best
    # site, so a floor would then ask "was the strongest site accessible?" instead of "has
    # this pair an accessible site at all?" — up to 3.5x fewer surviving pairs on four
    # genomes. NaN Pu fails the comparison, so unknown accessibility is never admitted.
    for column, floor in floors:
        if floor is not None:
            chunk = chunk[chunk[column] >= float(floor)]

    # Annotation runs on what survived the gate and floors, not on the majority that did
    # not. Position is a 0-1 fraction for cross-arm comparability with tidy_rnahybrid, which anchors the same way.
    return chunk.assign(
        Gene_length=lambda d: d["Gene"].map(gene_lengths).astype(int),
        Position=lambda d: d["Start1"].astype(float) / d["Gene_length"],
    )


#----- The per-pair cap, hoisted into the read loop -----#
def _precap_chunk(chunk, max_suboptimal_hits):
    if max_suboptimal_hits is None:
        return chunk

    ranked = chunk.sort_values("E_hybrid", kind="stable")
    return ranked.groupby(CAP_KEYS, sort=False).head(int(max_suboptimal_hits))


#----- Reads IntaRNA's CSV in batches, keeping only gated/floored rows, then ranks and caps on E_hybrid and emits the CSV -----#
def tidy_intarna(input_path, target_fasta_path, output_path,
                 max_hybrid_energy=None, max_suboptimal_hits=None,
                 min_target_unpaired_probability=None,
                 min_query_unpaired_probability=None):
    floors = (("Pu1", min_target_unpaired_probability),
              ("Pu2", min_query_unpaired_probability))
    # E/E_hybrid typed by the C parser rather than cast afterwards: an astype() would run on
    # the whole batch, most of which the gate and cap below are about to discard.
    read_options = dict(sep=";", dtype={"id1": str, "id2": str, "E": float, "E_hybrid": float},
                        na_values=INTARNA_NA_VALUES)

    header = pd.read_csv(input_path, nrows=0, **read_options)
    _require_intarna_columns(header.columns, floors)

    gene_lengths = _parse_gene_lengths(target_fasta_path)

    # Batched, and capped per batch: how much of a batch dies in the filters alone depends entirely on the config
    prepared = [
        _precap_chunk(_prepare_chunk(chunk, gene_lengths, max_hybrid_energy, floors),
                      max_suboptimal_hits)
        for chunk in pd.read_csv(input_path, chunksize=READ_CHUNK_ROWS, **read_options)
    ]
    df = pd.concat(prepared, ignore_index=True)

    # Rank by the gated quantity, then cap. Capping first could discard a qualifying row
    # in favour of a better-total-E one that fails the gate.
    df = df.sort_values("E_hybrid", kind="stable")

    if max_suboptimal_hits is not None:
        df = df.groupby(CAP_KEYS, sort=False).head(int(max_suboptimal_hits))

    # After the cap, not before: the correction ranks a pair's genes against each other, so a
    # gene should weigh as close to once as the cap allows. Column only, never a gate — see
    # _length.py.
    df = add_length_corrected(df, "E_hybrid")

    present_cols = [c for c in OUTPUT_COLUMNS if c in df.columns]
    df = df[present_cols]

    df.to_csv(ensure_parent(output_path), index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    tidy_intarna(
        input_path=snakemake.input.csv,
        target_fasta_path=snakemake.input.target,
        output_path=snakemake.output.tidy,
        max_hybrid_energy=snakemake.params.max_hybrid_energy,
        max_suboptimal_hits=snakemake.params.max_suboptimal_hits,
        min_target_unpaired_probability=snakemake.params.min_target_unpaired_probability,
        min_query_unpaired_probability=snakemake.params.min_query_unpaired_probability,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
