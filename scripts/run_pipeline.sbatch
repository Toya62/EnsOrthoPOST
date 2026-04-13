#!/bin/bash
#SBATCH --job-name=ensorthopost_test
#SBATCH --output=/work/%u/GLUH/EnsOrthoPOST/logs/ensorthopost_%j.log
#SBATCH --error=/work/%u/GLUH/EnsOrthoPOST/logs/ensorthopost_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

#SBATCH --mail-type=END,FAIL


mkdir -p /work/$USER/GLUH/EnsOrthoPOST/logs

set -euo pipefail

echo "========================================"
echo "EnsOrthoPOST Pipeline"
echo "Start: $(date)"
echo "========================================"
echo ""

cd /work/$USER/GLUH/EnsOrthoPOST

GENE_FAMILY_FILE="/work/$USER/GLUH/EnsOrthoPOST/gene_family_ids.txt"
OUTPUT_DIR="/work/$USER/GLUH/EnsOrthoPOST/output"
SPECIES_FILE="/work/$USER/GLUH/EnsOrthoPOST/species_list.txt"

if [[ ! -r "$GENE_FAMILY_FILE" ]]; then
    echo "ERROR: $GENE_FAMILY_FILE not found"
    exit 1
fi

# Parse train_group_ids.txt:
# - one numeric group ID per non-empty line
# - strips comments (#)
GENE_FAMILIES=()
while IFS= read -r line || [[ -n "$line" ]]; do
    id_clean="${line%%#*}"
    id_clean="${id_clean//[[:space:]]/}"
    if [[ -z "$id_clean" ]]; then
        continue
    fi
    if [[ "$id_clean" =~ ^[0-9]+$ ]]; then
        GENE_FAMILIES+=("$id_clean")
    else
        echo "ERROR: Invalid ID '$id_clean' in $GENE_FAMILY_FILE"
        exit 1
    fi
done < "$GENE_FAMILY_FILE"

TOTAL=${#GENE_FAMILIES[@]}
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No valid gene families found in $GENE_FAMILY_FILE"
    exit 1
fi

echo "Found $TOTAL test gene families to process"
echo "Output directory: $OUTPUT_DIR"
echo "Species file: $SPECIES_FILE"
echo "Human sequences in FASTA/GTF: excluded"
echo ""


for i in "${!GENE_FAMILIES[@]}"; do
    GROUP_ID="${GENE_FAMILIES[$i]}"
    CURRENT=$((i + 1))

    echo "========================================"
    echo "Processing $CURRENT/$TOTAL: Gene Family $GROUP_ID"
    echo "Start: $(date)"
    echo "========================================"

    set +e
    python3 main_orchestrator.py \
        --group-id "$GROUP_ID" \
        --workers 4 \
        --throttle 0.2 \
        --species-file "$SPECIES_FILE" \
        --output-dir "$OUTPUT_DIR" \
        --deterministic
    EXIT_CODE=$?
    set -e

    if [ $EXIT_CODE -eq 0 ]; then
        echo "OK: Family $GROUP_ID done: $(date)"
    elif [ $EXIT_CODE -eq 2 ]; then
        echo "WARNING: Family $GROUP_ID SKIPPED (no orthologs): $(date)"
    else
        echo "ERROR: Family $GROUP_ID FAILED (exit $EXIT_CODE): $(date)"
    fi
    echo "========================================"
    echo ""
done

echo "========================================"
echo "All $TOTAL families processed"
echo "End: $(date)"
echo "========================================"