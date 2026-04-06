#!/bin/bash
#PBS -N repeatmasker_chunk
#PBS -l select=1:ncpus=16:mem=192gb:scratch_local=200gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -v chunk

# run example:
# cd /storage/brno12-cerit/home/kratka/brno1/LTR-TE_dynamics/
# qsub workflow/scripts/run_repeatmasker_chunk.sh -- data/chromosomes_chunked/Triticum_monococcum/chr1/chunk_299700000.fasta

# used as a manual fallback for testing and failed snakemake jobs

# Set this to the root of your project
PROJECT_DIR="/storage/brno12-cerit/home/kratka/brno1/LTR-TE_dynamics"

set -euo pipefail

# === Environment Setup ===
module add repeatmasker
export TMPDIR="$SCRATCHDIR"
export OMP_NUM_THREADS="$PBS_NUM_PPN"

# === Paths ===
cd "$PROJECT_DIR"

CHUNK_FASTA="$chunk"

if [[ ! -f "$CHUNK_FASTA" ]]; then
    echo "❌ Error: input file '$CHUNK_FASTA' not found."
    exit 1
fi

# Extract relative path, chunk ID, and sample
REL_PATH="${CHUNK_FASTA#data/chromosomes_chunked/}"             # Triticum_monococcum/chr1/chunk_299700000.fasta
CHUNK_ID="${REL_PATH%.fasta}"                                   # Triticum_monococcum/chr1/chunk_299700000
SAMPLE="${REL_PATH%%/*}"                                        # Triticum_monococcum

# Derived paths
LIBRARY="$PROJECT_DIR/data/ltr_library/${SAMPLE}/TE_DLplus.fasta"
OUTPUT_DIR="$PROJECT_DIR/data/repeatmasker_chunked/$(dirname "$CHUNK_ID")"
mkdir -p "$OUTPUT_DIR"

# Prepare scratch
SCRATCH_FASTA="$SCRATCHDIR/$(basename "$CHUNK_FASTA")" 
cp "$CHUNK_FASTA" "$SCRATCH_FASTA"

# === Run RepeatMasker in scratch ===
cd "$SCRATCHDIR"
echo "🚀 Running RepeatMasker on $SCRATCH_FASTA"

RepeatMasker -pa "$OMP_NUM_THREADS" \
             -lib "$LIBRARY" \
             -gff \
             -dir . \
             "$SCRATCH_FASTA" \
             -no_is -nolow

# === Copy results back ===
CHUNK_BASE=$(basename "$SCRATCH_FASTA")

echo "📦 Copying all RepeatMasker output files back to $OUTPUT_DIR and cleaning scratch"
cp "${CHUNK_BASE}".* "$OUTPUT_DIR/" && rm -rf "$SCRATCHDIR"/*


echo "✅ Done for $CHUNK_FASTA"