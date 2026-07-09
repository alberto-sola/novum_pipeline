#!/usr/bin/env bash
#
# setup-macos.sh — one-time setup so the pipeline runs on Apple-Silicon Macs.
#
# RNAhybrid has no osx-arm64 conda build, so on Apple Silicon we resolve the whole
# conda toolchain as osx-64 (Intel) and run it under Rosetta 2. Run this ONCE, with
# your snakemake conda env active. No-op on Linux and Intel Macs.
#
set -euo pipefail

os="$(uname -s)"
machine="$(uname -m)"

# --- 1. Platform guard ---
if [ "$os" != "Darwin" ]; then
    echo "setup-macos.sh: not macOS (detected '$os') — nothing to do."
    exit 0
fi
if [ "$machine" != "arm64" ]; then
    echo "setup-macos.sh: Intel Mac (detected '$machine') — no setup needed; conda resolves osx-64 natively."
    exit 0
fi

echo "setup-macos.sh: Apple Silicon detected — configuring osx-64 + Rosetta 2."

# --- 2. Require an active, non-base conda env ---
activate_hint="Activate your snakemake env first, e.g.:  conda activate snakemake-modern"
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: no active conda environment detected (CONDA_PREFIX is unset)." >&2
    echo "$activate_hint" >&2
    exit 1
fi
if [ "${CONDA_DEFAULT_ENV:-}" = "base" ]; then
    echo "ERROR: the conda 'base' environment is active." >&2
    echo "$activate_hint" >&2
    exit 1
fi

# --- 3. Ensure Rosetta 2 ---
rosetta_install=(softwareupdate --install-rosetta --agree-to-license)
if arch -x86_64 /usr/bin/true >/dev/null 2>&1; then
    echo "Rosetta 2: already installed."
else
    echo "Rosetta 2: not found — installing (may prompt for your password)..."
    if ! "${rosetta_install[@]}"; then
        echo "ERROR: Rosetta 2 installation failed." >&2
        echo "Install it manually with:  ${rosetta_install[*]}" >&2
        exit 1
    fi
fi

# --- 4. Persist the osx-64 subdir in THIS env ---
conda config --env --set subdir osx-64
echo "Pinned conda subdir to osx-64 for env: ${CONDA_DEFAULT_ENV:-$CONDA_PREFIX}"

# --- 5. Confirm + next step ---
cat <<'EOF'

Done — this env now resolves all conda packages as osx-64 (Intel) under Rosetta 2.

Next:
  snakemake --use-conda --cores 8      # or launch the UI and click Run

Notes:
  * The FIRST run downloads/builds Intel packages and is slower; later runs are normal.
  * One-time step per env — re-run setup-macos.sh only if you recreate the env.
EOF
