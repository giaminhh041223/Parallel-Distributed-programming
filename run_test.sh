#!/usr/bin/env bash

# ============================================================================
#  run_test.sh — Automated Compilation, Validation & Benchmarking Script
#
#  Usage:
#    ./run_test.sh [hostfile]              # Default block_size=128
#    ./run_test.sh [hostfile] [block_size] # Custom block size
#
#  Features:
#    1. Compiles seq_align (g++) and mpi_align (mpicxx/mpic++/mpicc)
#    2. Correctness verification at 5000×5000
#    3. Multi-size, multi-process scalability benchmark
#    4. CSV output for direct import into Python plotting scripts
#    5. Per-rank timing breakdown (load balance data)
# ============================================================================

set -e  # Exit on first error

# ── Color helpers for terminal output ──
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'  # No Color

# ── Parse arguments ──
HOSTFILE_ARG=""
BLOCK_SIZE=128

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 [hostfile] [block_size]"
    echo "  hostfile  : (Optional) MPI hostfile for distributed execution"
    echo "  block_size: (Optional) Block size B for MPI version (default: 128)"
    exit 0
fi

if [ -n "$1" ]; then
    if [ ! -f "$1" ]; then
        echo -e "${RED}Error: Hostfile '$1' not found!${NC}"
        exit 1
    fi
    HOSTFILE_ARG="--hostfile $1"
    echo -e "${CYAN}Using MPI Hostfile: $1${NC}"
fi

if [ -n "$2" ]; then
    BLOCK_SIZE=$2
fi

echo -e "${BOLD}Block Size (B): ${BLOCK_SIZE}${NC}"
echo ""

# ── MPI run flags ──
# --oversubscribe: allows more processes than physical cores (for testing)
# --allow-run-as-root: needed on some cluster setups
MPI_FLAGS="--oversubscribe"
if [ "$(id -u)" -eq 0 ]; then
    MPI_FLAGS="$MPI_FLAGS --allow-run-as-root"
fi

# ===================================================================
#  1. COMPILATION
# ===================================================================
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  PHASE 1: COMPILATION${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"

echo -e "${CYAN}Compiling sequential baseline (g++ -O3)...${NC}"
g++ -O3 -o seq_align src/seq_align.cpp

echo -e "${CYAN}Compiling MPI parallel version (mpicxx -O3)...${NC}"
if command -v mpicxx >/dev/null 2>&1; then
    mpicxx -O3 -o mpi_align src/mpi_align.cpp
elif command -v mpic++ >/dev/null 2>&1; then
    mpic++ -O3 -o mpi_align src/mpi_align.cpp
else
    mpicc -O3 -o mpi_align src/mpi_align.cpp -lstdc++
fi
echo -e "${GREEN}Compilation complete!${NC}"
echo ""

# ===================================================================
#  2. CORRECTNESS VERIFICATION (5000 × 5000)
# ===================================================================
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  PHASE 2: CORRECTNESS VERIFICATION (5000×5000)${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"

echo -e "${CYAN}Running Sequential Baseline (5000x5000)...${NC}"
SEQ_OUTPUT=$(./seq_align 5000 5000 42)
SEQ_SCORE=$(echo "$SEQ_OUTPUT" | grep "Max Alignment Score:" | awk '{print $NF}')
echo "  Sequential Score: $SEQ_SCORE"

echo -e "${CYAN}Running MPI (np=3, B=$BLOCK_SIZE, 5000x5000)...${NC}"
MPI_OUTPUT=$(mpirun -np 3 $MPI_FLAGS $HOSTFILE_ARG ./mpi_align 5000 5000 42 $BLOCK_SIZE)
MPI_SCORE=$(echo "$MPI_OUTPUT" | grep "Global Max Alignment Score:" | awk '{print $NF}')
echo "  MPI Parallel Score: $MPI_SCORE"

if [ "$SEQ_SCORE" != "$MPI_SCORE" ]; then
    echo -e "\n${RED}${BOLD}  ❌ VERIFICATION FAILED for 5000x5000!${NC}"
    echo -e "  Sequential Score: ${SEQ_SCORE}"
    echo -e "  MPI Score:        ${MPI_SCORE}"
    exit 1
fi

echo -e "${CYAN}Running Sequential Baseline (7777x7777)...${NC}"
SEQ_OUTPUT2=$(./seq_align 7777 7777 42)
SEQ_SCORE2=$(echo "$SEQ_OUTPUT2" | grep "Max Alignment Score:" | awk '{print $NF}')
echo "  Sequential Score: $SEQ_SCORE2"

echo -e "${CYAN}Running MPI (np=3, B=$BLOCK_SIZE, 7777x7777)...${NC}"
MPI_OUTPUT2=$(mpirun -np 3 $MPI_FLAGS $HOSTFILE_ARG ./mpi_align 7777 7777 42 $BLOCK_SIZE)
MPI_SCORE2=$(echo "$MPI_OUTPUT2" | grep "Global Max Alignment Score:" | awk '{print $NF}')
echo "  MPI Parallel Score: $MPI_SCORE2"

if [ "$SEQ_SCORE2" != "$MPI_SCORE2" ]; then
    echo -e "\n${RED}${BOLD}  ❌ VERIFICATION FAILED for 7777x7777!${NC}"
    echo -e "  Sequential Score: ${SEQ_SCORE2}"
    echo -e "  MPI Score:        ${MPI_SCORE2}"
    exit 1
fi

echo -e "\n${GREEN}${BOLD}  ✅ VERIFICATION PASSED: All scores match exactly!${NC}\n"

# ===================================================================
#  3. SCALABILITY BENCHMARKS
# ===================================================================
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  PHASE 3: SCALABILITY BENCHMARKS${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"

SIZES=(5000 10000 20000)
NP_LIST=(1 2 3 4 6 8 12)

# Create CSV output file
CSV_FILE="benchmark_results.csv"
echo "size,num_procs,time_ms,score" > "$CSV_FILE"

# Header for terminal table
printf "\n${BOLD}%-12s" "Size"
for NP in "${NP_LIST[@]}"; do
    printf "| np=%-8s" "$NP"
done
printf "| Sequential${NC}\n"
printf "%-12s" "────────────"
for NP in "${NP_LIST[@]}"; do
    printf "|────────────"
done
printf "|────────────\n"

for SIZE in "${SIZES[@]}"; do
    printf "%-12s" "${SIZE}×${SIZE}"

    for NP in "${NP_LIST[@]}"; do
        MPI_OUT=$(mpirun -np $NP $MPI_FLAGS $HOSTFILE_ARG ./mpi_align $SIZE $SIZE 42 $BLOCK_SIZE 2>/dev/null)
        MPI_TIME=$(echo "$MPI_OUT" | grep "Parallel Execution Time:" | awk '{print $(NF-1)}')
        MPI_SC=$(echo "$MPI_OUT" | grep "Global Max Alignment Score:" | awk '{print $NF}')
        echo "$SIZE,$NP,$MPI_TIME,$MPI_SC" >> "$CSV_FILE"
        printf "| %8s ms" "$MPI_TIME"
    done

    # Sequential
    SEQ_OUT=$(./seq_align $SIZE $SIZE 42)
    SEQ_TIME=$(echo "$SEQ_OUT" | grep "Execution Time:" | awk '{print $(NF-1)}')
    SEQ_SC=$(echo "$SEQ_OUT" | grep "Max Alignment Score:" | awk '{print $NF}')
    echo "$SIZE,0,$SEQ_TIME,$SEQ_SC" >> "$CSV_FILE"
    printf "| %8s ms\n" "$SEQ_TIME"
done

echo ""
echo -e "${GREEN}${BOLD}Benchmark results saved to: ${CSV_FILE}${NC}"

# ===================================================================
#  4. PER-RANK LOAD BALANCE (largest size, np=3)
# ===================================================================
echo ""
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  PHASE 4: LOAD BALANCE ANALYSIS (np=3)${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════${NC}"

LARGEST="${SIZES[${#SIZES[@]}-1]}"
echo -e "${CYAN}Running MPI (np=3, size=${LARGEST}×${LARGEST}, B=$BLOCK_SIZE)...${NC}"
mpirun -np 3 $MPI_FLAGS $HOSTFILE_ARG ./mpi_align $LARGEST $LARGEST 42 $BLOCK_SIZE

echo ""
echo -e "${GREEN}${BOLD}════════════════════════════════════════════${NC}"
echo -e "${GREEN}${BOLD}  ALL BENCHMARKS COMPLETE${NC}"
echo -e "${GREEN}${BOLD}  CSV: ${CSV_FILE}${NC}"
echo -e "${GREEN}${BOLD}════════════════════════════════════════════${NC}"
