#!/bin/bash
# run_demo_fasta.sh
# Run FASTA alignment demo on the MPI cluster.
# Usage: ./run_demo_fasta.sh [hostfile]

HOSTFILE=$1
BLOCK_SIZE=128

echo "=========================================================="
echo "          SMITH-WATERMAN FASTA LOCAL ALIGNMENT DEMO       "
echo "=========================================================="
echo "Dataset: SARS-CoV-2 (29903 bp) vs MERS-CoV (30119 bp)"
echo ""

echo "[STEP 1] Running Sequential Baseline..."
echo "----------------------------------------------------------"
./seq_align data/sars_cov_2.fasta data/mers_cov.fasta
echo ""

echo "[STEP 2] Running Parallel MPI Version..."
echo "----------------------------------------------------------"
if [ -z "$HOSTFILE" ]; then
    echo "No hostfile provided. Running locally on 3 ranks..."
    mpirun -np 3 ./mpi_align data/sars_cov_2.fasta data/mers_cov.fasta $BLOCK_SIZE
else
    echo "Using hostfile: $HOSTFILE"
    mpirun -np 12 --hostfile "$HOSTFILE" ./mpi_align data/sars_cov_2.fasta data/mers_cov.fasta $BLOCK_SIZE
fi
echo ""
echo "=========================================================="
echo "Demo completed."
echo "=========================================================="
