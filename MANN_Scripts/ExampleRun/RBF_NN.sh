#!/bin/bash
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --job-name=RBFNN

module load miniconda3
source activate tf_env
cd "$SLURM_SUBMIT_DIR"

## Change trial number to edit file naming
TRIAL=1


TASK_NAME=${SLURM_JOB_NAME:-RBFNN}
DATE=$(date +%F_%H-%M)
TAG=${TASK_NAME}_T${TRIAL}_${DATE}
[[ -n "$SLURM_ARRAY_TASK_ID" ]] && TAG=${TAG}_A${SLURM_ARRAY_TASK_ID}

OUTDIR="$SLURM_SUBMIT_DIR/$TAG"
mkdir -p "$OUTDIR"
export TRIAL TASK_NAME RUN_TAG="$TAG" OUTDIR

cd "$OUTDIR"

## script expects python script to be named 'RBF_NN.py'
## this can be changed if needed

python "$SLURM_SUBMIT_DIR/RBF_NN.py" > "$TAG.out" 2> "$TAG.err"

rm "slurm-${SLURM_JOB_ID}.out"



