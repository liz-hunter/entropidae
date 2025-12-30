#!/usr/bin/env bash
#SBATCH -J kraken_entropy
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --array=1-6
#SBATCH --output=%x_%A_%a.out
#SBATCH --time=24:00:00

module load python/3.12.5  
module load micromamba
micromamba activate ~/micromamba/envs/virid_env
micromamba shell reinit --shell bash

echo "START"
echo "JobID: ${SLURM_JOB_ID}  ArrayTaskID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
date

# define variables
WORKDIR="$MYPATH/entropy_test"
cd "$WORKDIR"

NODES_DMP="$MYPATH/entropy_stan/stan_full/taxonomy/nodes.dmp"
OUTDIR="$MYPATH/entropy_test"
PERREAD_OUTDIR="$OUTDIR/per_read"
mkdir -p "$OUTDIR" "$PERREAD_OUTDIR"

# split files into chunks with awk prior to running
CHUNKFILE="$WORKDIR/chunk_${SLURM_ARRAY_TASK_ID}.txt"
if [[ ! -f "$CHUNKFILE" ]]; then
  echo "ERROR: chunk file not found: $CHUNKFILE" >&2
  exit 2
fi

# Internal multiprocessing per array task (<= -c)
PY_WORKERS="${PY_WORKERS:-8}"

ALPHA_UP="${ALPHA_UP:-0.3}" #default
ALPHA_DOWN="${ALPHA_DOWN:-1.0}" #default
COMPRESSLEVEL="${COMPRESSLEVEL:-3}"

echo "Chunkfile: $CHUNKFILE"
echo "Python workers: $PY_WORKERS"
echo "Outdir: $PERREAD_OUTDIR"

python "$WORKDIR/weighted_entropy_batch.py" \
  --nodes-dmp "$NODES_DMP" \
  --jobs-tsv "$CHUNKFILE" \
  --outdir "$PERREAD_OUTDIR" \
  --summary-tsv "$OUTDIR/summary_chunk_${SLURM_ARRAY_TASK_ID}.tsv" \
  --jobs "$PY_WORKERS" \
  --alpha-up "$ALPHA_UP" \
  --alpha-down "$ALPHA_DOWN" \
  --gzip \
  --compresslevel "$COMPRESSLEVEL" \
  --skip-existing \
  --no-diagnostics

echo "END"
date

