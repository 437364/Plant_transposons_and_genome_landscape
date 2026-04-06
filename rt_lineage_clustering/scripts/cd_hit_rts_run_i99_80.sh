#!/bin/bash
#PBS -N cdhit_HordW_Pnot_Phal_Hvul_RTs_allTeFastas_PROT
#PBS -l select=1:ncpus=48:mem=100gb:scratch_local=20gb
#PBS -l walltime=24:00:00
#PBS -m abe

# === User variables ===
DATADIR=/storage/brno12-cerit/home/pj_narutro/HordeumSegments/HordW_Hvul_Phal_Pnot_RTs/merged_fastas
OUTDIR="$DATADIR"/CD_HIT_i99_80_PROT_PBS
mkdir -p "$OUTDIR"

# === Job info ===
echo "$PBS_JOBID is running on node $(hostname -f) in scratch directory $SCRATCHDIR" >> "$DATADIR/jobs_info.txt"

# === Load CD-HIT ===
module load cdhit/4.8.1

# === Check scratch ===
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# === Copy all FASTAs to scratch ===
cp "$DATADIR"/*_all_RTs.fa "$SCRATCHDIR"/ || { echo >&2 "Error copying input FASTAs!"; exit 2; }
cd "$SCRATCHDIR" || exit 1

# === Identity thresholds (protein mode) ===
# Note: valid range for cd-hit (protein) is typically 0.4–1.0
thresholds=(0.99 0.98 0.97 0.96 0.95 0.93 0.90 0.85 0.80)

# === Loop over each FASTA file ===
for fasta in *_all_RTs.fa; do
    pref=$(basename "$fasta" .fa)
    echo "=== Processing file: $pref ==="

    # Loop over all thresholds
    for c in "${thresholds[@]}"; do
        # Set valid word size (-n) for protein mode
        # Reference: cd-hit manual
        if (( $(echo "$c >= 0.7" | bc -l) )); then
            n=5
        elif (( $(echo "$c >= 0.6" | bc -l) )); then
            n=4
        else
            n=3
        fi

        echo "Running cd-hit (protein) on ${pref}.fa with -c $c and -n $n"
        cd-hit -T 48 -M 100000 -i "$fasta" -o "${pref}_c${c//./}.fa" -c "$c" -n "$n"

        # Check success and copy results
        if [[ $? -eq 0 ]]; then
            cp "${pref}_c${c//./}.fa"* "$OUTDIR"/ || { echo >&2 "Copying results failed for $pref (c=$c)!"; }
        else
            echo "❌ cd-hit failed for $pref (c=$c)" >&2
        fi
    done
done

# === Clean up scratch ===
clean_scratch
