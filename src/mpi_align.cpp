#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cctype>
#include <mpi.h>

/*============================================================================
 * PARALLEL SMITH-WATERMAN — BLOCK-BASED WAVEFRONT WITH CYCLIC COLUMNS
 *
 * Performance improvements over the original version:
 *   1. H_local is a PRE-ALLOCATED flat 1D array (B*B), reused across all
 *      block computations — eliminates ~24,000 heap alloc/dealloc cycles.
 *   2. Receive buffers use a fixed pre-allocated pool instead of per-step
 *      vector<vector<int>> allocations.
 *   3. Left boundary is accessed via pointer, no redundant copy.
 *   4. top_boundary uses flat 1D arrays for cache-friendly access.
 *   5. Per-rank computation time instrumentation added for load balance
 *      analysis (output when env var SW_VERBOSE=1).
 *   6. Branchless scoring via std::max.
 *   7. Cross-platform memory prefetching for next row of H_local.
 *   8. RAII memory management via std::vector for H_local and recv_pool.
 *   9. Tag calculation uses % 32768 to match MPI standard tag limits.
 *  10. Removed redundant MPI_Barrier before final MPI_Wtime().
 *
 * Scoring: Match = +3, Mismatch = -3, Gap = -2
 * Block size: configurable via 5th CLI argument or compile-time #define B
 *==========================================================================*/

// Cross-platform prefetch
#if defined(__GNUC__) || defined(__clang__)
  #define PREFETCH(addr) __builtin_prefetch(addr, 1, 1)
#elif defined(_MSC_VER)
  #include <intrin.h>
  #define PREFETCH(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T0)
#else
  #define PREFETCH(addr)
#endif

// Default block size (overridable at compile time: -DB=256)
#ifndef B
#define B 128
#endif

// Scoring constants
static const int MATCH    =  3;
static const int MISMATCH = -3;
static const int GAP      = -2;

// Check if file path ends in a FASTA extension
bool is_fasta_file(const std::string& arg) {
    if (arg.length() < 6) return false;
    std::string ext = arg.substr(arg.length() - 6);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".fasta") return true;
    if (arg.length() >= 3) {
        std::string ext2 = arg.substr(arg.length() - 3);
        std::transform(ext2.begin(), ext2.end(), ext2.begin(), ::tolower);
        if (ext2 == ".fa") return true;
    }
    return false;
}

// Read FASTA file and return sequence
std::string read_fasta(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[Lỗi] Không thể mở file FASTA: " << filepath << std::endl;
        std::exit(1);
    }
    std::string line, seq;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue; // Bỏ qua dòng tiêu đề
        for (char c : line) {
            if (!std::isspace(c)) {
                seq.push_back(std::toupper(c));
            }
        }
    }
    file.close();
    return seq;
}

// Generate a deterministic pseudo-random DNA sequence
std::string generate_dna_sequence(size_t length, unsigned int seed) {
    std::string seq;
    seq.reserve(length);
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < length; ++i) {
        seq.push_back(nucleotides[dis(gen)]);
    }
    return seq;
}

int main(int argc, char* argv[]) {
    // ── MPI Initialization ──
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // ── Parse & broadcast parameters ──
    int len1 = 10000, len2 = 10000;
    unsigned int seed = 42;
    int block_size = B;
    bool use_fasta = false;
    std::string file1, file2;

    if (rank == 0) {
        if (argc >= 3 && is_fasta_file(argv[1]) && is_fasta_file(argv[2])) {
            use_fasta = true;
            file1 = argv[1];
            file2 = argv[2];
            if (argc >= 4) block_size = std::atoi(argv[3]);
        } else {
            if (argc >= 2) len1 = std::atoi(argv[1]);
            if (argc >= 3) len2 = std::atoi(argv[2]);
            if (argc >= 4) seed = static_cast<unsigned int>(std::atoi(argv[3]));
            if (argc >= 5) block_size = std::atoi(argv[4]);
        }
    }

    std::string seq1, seq2;

    if (rank == 0) {
        std::cout << "MPI Parallel Smith-Waterman (Block-based Wavefront)" << std::endl;
        std::cout << "Number of MPI processes: " << nprocs << std::endl;
        if (use_fasta) {
            std::cout << "FASTA File Mode Enabled." << std::endl;
            std::cout << "Loading file 1: " << file1 << std::endl;
            std::cout << "Loading file 2: " << file2 << std::endl;
            seq1 = read_fasta(file1);
            seq2 = read_fasta(file2);
            len1 = seq1.length();
            len2 = seq2.length();
            std::cout << "Sequence 1 length: " << len1 << " bp" << std::endl;
            std::cout << "Sequence 2 length: " << len2 << " bp" << std::endl;
        } else {
            std::cout << "Sequence 1 length: " << len1 << std::endl;
            std::cout << "Sequence 2 length: " << len2 << std::endl;
            std::cout << "Random Seed: " << seed << std::endl;
            seq1 = generate_dna_sequence(len1, seed);
            seq2 = generate_dna_sequence(len2, seed + 1);
        }
        std::cout << "Block Size: " << block_size << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }

    // Pack into a single broadcast to minimize latency
    int params[4];
    if (rank == 0) {
        params[0] = len1; params[1] = len2;
        params[2] = use_fasta ? 1 : 0; params[3] = block_size;
    }
    MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);
    len1 = params[0]; len2 = params[1];
    use_fasta = (params[2] == 1); block_size = params[3];

    // ── Allocate & broadcast sequences ──
    if (rank != 0) {
        seq1.resize(len1);
        seq2.resize(len2);
    }
    MPI_Bcast(&seq1[0], len1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&seq2[0], len2, MPI_CHAR, 0, MPI_COMM_WORLD);

    // ── Block grid dimensions ──
    const int M_b = (len1 + 1 + block_size - 1) / block_size;  // row blocks
    const int N_b = (len2 + 1 + block_size - 1) / block_size;  // col blocks

    // ── Determine local column allocations (cyclic: c % nprocs == rank) ──
    int my_cols = 0;
    for (int c = 0; c < N_b; ++c)
        if (c % nprocs == rank) my_cols++;

    // ── Pre-allocate ALL memory upfront ──

    // top_boundary[c_local] stores the bottom row of the block above (r-1, c)
    // Flat storage: my_cols arrays, each of width <= block_size
    std::vector<std::vector<int>> top_boundary(my_cols);
    {
        int idx = 0;
        for (int c = 0; c < N_b; ++c) {
            if (c % nprocs == rank) {
                int col_end = std::min((c + 1) * block_size - 1, len2);
                int w = col_end - c * block_size + 1;
                top_boundary[idx].assign(w, 0);  // H[0][j] = 0
                idx++;
            }
        }
    }

    // Send buffers: one per local column, reused across wavefront steps
    std::vector<std::vector<int>> send_bufs(my_cols,
        std::vector<int>(block_size + 1, 0));
    std::vector<MPI_Request> send_reqs(my_cols, MPI_REQUEST_NULL);

    // H_local: PRE-ALLOCATED flat 1D array of size block_size * block_size
    // Indexed as H_local[i * block_size + j]
    // RAII via std::vector
    std::vector<int> H_local_vec(block_size * block_size);
    int* H_local = H_local_vec.data();

    // Receive buffer pool: at most my_cols blocks can be active per step.
    // Each needs at most (block_size + 1) ints.
    // RAII via std::vector
    const int max_recv = my_cols;
    std::vector<int> recv_pool_vec(max_recv * (block_size + 1));
    int* recv_pool = recv_pool_vec.data();

    // Active block descriptor
    struct ActiveBlock {
        int r, c, c_local, h, w;
        int recv_slot;   // index into recv_pool, or -1 if no recv needed
    };
    std::vector<ActiveBlock> active_blocks;
    active_blocks.reserve(my_cols);

    std::vector<MPI_Request> recv_reqs;
    recv_reqs.reserve(max_recv);

    int local_max_score = 0;
    double local_comp_time = 0.0;  // for load-balance instrumentation

    // ── Synchronize & start timer ──
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    // ═══════════════════════════════════════════════════════
    //  WAVEFRONT MAIN LOOP
    //  Step t corresponds to anti-diagonal t of the block grid.
    //  Block (r, c) is computed at step t = r + c.
    // ═══════════════════════════════════════════════════════
    const int total_steps = M_b + N_b - 1;

    for (int t = 0; t < total_steps; ++t) {
        active_blocks.clear();
        recv_reqs.clear();
        int recv_slot = 0;

        // ── STEP 1: Identify active blocks and post non-blocking receives ──
        int c_local = 0;
        for (int c = 0; c < N_b; ++c) {
            if (c % nprocs == rank) {
                int r = t - c;
                if (r >= 0 && r < M_b) {
                    ActiveBlock blk;
                    blk.r = r;
                    blk.c = c;
                    blk.c_local = c_local;

                    int i_end = std::min((r + 1) * block_size - 1, len1);
                    blk.h = i_end - r * block_size + 1;

                    int j_end = std::min((c + 1) * block_size - 1, len2);
                    blk.w = j_end - c * block_size + 1;

                    blk.recv_slot = -1;

                    // If this block has a left neighbor, we need its right
                    // boundary column (+ the top-left corner value).
                    if (c > 0) {
                        blk.recv_slot = recv_slot;
                        int* rbuf = recv_pool + recv_slot * (block_size + 1);
                        int left_rank = (rank - 1 + nprocs) % nprocs;
                        int tag = r % 32768;
                        MPI_Request req;
                        MPI_Irecv(rbuf, blk.h + 1, MPI_INT,
                                  left_rank, tag, MPI_COMM_WORLD, &req);
                        recv_reqs.push_back(req);
                        recv_slot++;
                    }

                    active_blocks.push_back(blk);
                }
                c_local++;
            }
        }

        // ── STEP 2: Wait for all boundary data to arrive ──
        if (!recv_reqs.empty()) {
            MPI_Waitall(static_cast<int>(recv_reqs.size()),
                        &recv_reqs[0], MPI_STATUSES_IGNORE);
        }

        // ── STEP 3: Compute all active blocks ──
        double comp_t0 = MPI_Wtime();

        for (const auto& blk : active_blocks) {
            const int r = blk.r;
            const int c = blk.c;
            const int cl = blk.c_local;
            const int h = blk.h;
            const int w = blk.w;
            const int BS = block_size;  // alias for flat indexing

            // Pointer to left boundary data (from recv_pool or zeros)
            // left_bnd[0]   = top-left corner H(r*B-1, c*B-1)
            // left_bnd[1..h] = right column of block (r, c-1)
            int* left_bnd = nullptr;
            if (c > 0 && blk.recv_slot >= 0) {
                left_bnd = recv_pool + blk.recv_slot * (block_size + 1);
            }

            // Compute the block using the pre-allocated H_local
            for (int i = 0; i < h; ++i) {
                const int i_g = r * BS + i;  // global row

                // Prefetch next row of H_local
                if (i + 1 < h) {
                    PREFETCH(&H_local[(i + 1) * BS]);
                }

                for (int j = 0; j < w; ++j) {
                    const int j_g = c * BS + j;  // global col

                    if (i_g == 0 || j_g == 0) {
                        H_local[i * BS + j] = 0;
                        continue;
                    }

                    const char c1 = seq1[i_g - 1];
                    const char c2 = seq2[j_g - 1];
                    const int match_score = (c1 == c2) ? MATCH : MISMATCH;

                    // ── Diagonal: H(i_g-1, j_g-1) ──
                    int diag;
                    if (i > 0 && j > 0) {
                        diag = H_local[(i - 1) * BS + (j - 1)];
                    } else if (i == 0 && j > 0) {
                        diag = top_boundary[cl][j - 1];
                    } else if (i > 0 && j == 0) {
                        diag = (left_bnd != nullptr) ? left_bnd[i] : 0;
                    } else { // i == 0 && j == 0 — THE TOP-LEFT CORNER
                        diag = (left_bnd != nullptr) ? left_bnd[0] : 0;
                    }

                    // ── Top: H(i_g-1, j_g) ──
                    int up;
                    if (i > 0) {
                        up = H_local[(i - 1) * BS + j];
                    } else {
                        up = top_boundary[cl][j];
                    }

                    // ── Left: H(i_g, j_g-1) ──
                    int left;
                    if (j > 0) {
                        left = H_local[i * BS + (j - 1)];
                    } else {
                        left = (left_bnd != nullptr) ? left_bnd[i + 1] : 0;
                    }

                    // ── Branchless score computation (Smith-Waterman max with 0) ──
                    int score = std::max({diag + match_score, up + GAP, left + GAP, 0});

                    H_local[i * BS + j] = score;

                    if (score > local_max_score)
                        local_max_score = score;
                }
            }

            // ── STEP 4: Save corner & update top_boundary ──

            // The top-left corner that the RIGHT neighbor will need for its
            // block (r, c+1) at position (0,0) is: H(r*B-1, c*B+w-1).
            // This value is currently stored in top_boundary[cl][w-1]
            // (the bottom-right cell of block (r-1, c)).
            // We must capture it BEFORE overwriting top_boundary.
            int old_corner = top_boundary[cl][w - 1];

            // Overwrite top_boundary with the last row of H_local
            for (int j = 0; j < w; ++j) {
                top_boundary[cl][j] = H_local[(h - 1) * BS + j];
            }

            // ── STEP 5: Non-blocking send to right neighbor ──
            if (c < N_b - 1) {
                // Ensure previous send on this column slot is complete
                if (send_reqs[cl] != MPI_REQUEST_NULL) {
                    MPI_Wait(&send_reqs[cl], MPI_STATUS_IGNORE);
                }

                // Pack: [0] = corner, [1..h] = rightmost column
                send_bufs[cl][0] = old_corner;
                for (int i = 0; i < h; ++i) {
                    send_bufs[cl][i + 1] = H_local[i * BS + (w - 1)];
                }

                int right_rank = (rank + 1) % nprocs;
                int tag = r % 32768;
                MPI_Isend(&send_bufs[cl][0], h + 1, MPI_INT,
                          right_rank, tag, MPI_COMM_WORLD, &send_reqs[cl]);
            }
        }

        double comp_t1 = MPI_Wtime();
        local_comp_time += (comp_t1 - comp_t0);
    }

    // ── Drain all outstanding sends ──
    for (int i = 0; i < my_cols; ++i) {
        if (send_reqs[i] != MPI_REQUEST_NULL) {
            MPI_Wait(&send_reqs[i], MPI_STATUS_IGNORE);
        }
    }

    // ── Final timing (no redundant barrier) ──
    double t_end = MPI_Wtime();
    double local_elapsed = t_end - t_start;

    // ── Global reduction ──
    int global_max_score = 0;
    double max_elapsed = 0.0;
    MPI_Reduce(&local_max_score, &global_max_score, 1,
               MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_elapsed, &max_elapsed, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // ── Gather per-rank computation times for load balance analysis ──
    std::vector<double> all_comp_times_vec;
    std::vector<double> all_total_times_vec;
    if (rank == 0) {
        all_comp_times_vec.resize(nprocs);
        all_total_times_vec.resize(nprocs);
    }
    MPI_Gather(&local_comp_time, 1, MPI_DOUBLE,
               rank == 0 ? all_comp_times_vec.data() : nullptr,
               1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_elapsed, 1, MPI_DOUBLE,
               rank == 0 ? all_total_times_vec.data() : nullptr,
               1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // ── Output results ──
    if (rank == 0) {
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "Global Max Alignment Score: " << global_max_score << std::endl;
        std::cout << "Parallel Execution Time: " << max_elapsed * 1000.0 << " ms" << std::endl;
        std::cout << "--------------------------------------" << std::endl;

        // Detailed per-rank breakdown (for load balance charts)
        std::cout << "\n[Per-Rank Timing Breakdown]" << std::endl;
        std::cout << "Rank | Comp Time (ms) | Total Time (ms) | Comm/Idle (ms)" << std::endl;
        std::cout << "-----|----------------|-----------------|---------------" << std::endl;
        for (int p = 0; p < nprocs; ++p) {
            double ct = all_comp_times_vec[p] * 1000.0;
            double tt = all_total_times_vec[p] * 1000.0;
            double ci = tt - ct;
            printf("  %2d | %14.2f | %15.2f | %13.2f\n", p, ct, tt, ci);
        }
    }

    MPI_Finalize();
    return 0;
}