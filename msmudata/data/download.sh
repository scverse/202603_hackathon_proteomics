#!/bin/bash
set -euo pipefail

MINIMAL_URL="XXX"
ALBRECHT2025_URL="https://datashare.biochem.mpg.de/s/Cn2mfDbszYNS98e/download"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] [DATASET...]

Download proteomics datasets for msmudata.

Datasets:
  minimal         Small test dataset
  albrecht2025    Albrecht et al. 2025 study

If no dataset is specified, all datasets are downloaded.

Options:
  -o, --output DIR   Output directory (default: script directory)
  -h, --help         Show this help message
EOF
}

download_dataset() {
    local name="$1"
    local url="$2"
    local outdir="$3"
    local outfile="${outdir}/${name}.zip"

    if [[ "$url" == "XXX" ]]; then
        echo "Error: Download URL for '${name}' is not configured yet." >&2
        return 1
    fi

    echo "Downloading ${name}..."
    wget -O "$outfile" "$url"
    echo "Saved to ${outfile}"
}

outdir="$SCRIPT_DIR"
datasets=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)
            outdir="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        minimal|albrecht2025)
            datasets+=("$1")
            shift
            ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# Default to all datasets if none specified
if [[ ${#datasets[@]} -eq 0 ]]; then
    datasets=(minimal albrecht2025)
fi

mkdir -p "$outdir"

for ds in "${datasets[@]}"; do
    case "$ds" in
        minimal)      download_dataset minimal "$MINIMAL_URL" "$outdir" ;;
        albrecht2025) download_dataset albrecht2025 "$ALBRECHT2025_URL" "$outdir" ;;
    esac
done
