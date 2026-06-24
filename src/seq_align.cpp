#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cctype>

/*============================================================================
 * SEQUENTIAL SMITH-WATERMAN LOCAL ALIGNMENT — OPTIMIZED BASELINE
 *
 * Optimizations applied:
 *   1. Flat C-style array for zero-overhead indexing and contiguous memory.
 *   2. Single-row rolling buffer with explicit "diag" variable,
 *      eliminating the need for a second row array.
 *   3. Branchless std::max for score computation, enabling compiler to emit
 *      conditional moves (cmov) instead of branches.
 *   4. Loop-hoisted character fetch to reduce redundant loads.
 *   5. Compile with -O3 for auto-vectorization.
 *
 * Scoring: Match = +3, Mismatch = -3, Gap = -2
 *==========================================================================*/

static constexpr int MATCH    =  3;
static constexpr int MISMATCH = -3;
static constexpr int GAP      = -2;

// Generate a deterministic pseudo-random DNA sequence
std::string generate_dna_sequence(size_t length, unsigned int seed) {
    static constexpr char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string seq;
    seq.reserve(length);
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);
    for (size_t i = 0; i < length; ++i) {
        seq.push_back(nucleotides[dis(gen)]);
    }
    return seq;
}

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

/*--------------------------------------------------------------------------
 * Smith-Waterman using a SINGLE ROW + diagonal accumulator.
 * Space complexity: O(len2) instead of O(len1 * len2).
 *
 * prev[j]  holds H(i-1, j) at the start of row i.
 * diag     holds H(i-1, j-1) carried across the inner loop.
 *------------------------------------------------------------------------*/
int smith_waterman_sequential(const std::string& seq1, const std::string& seq2) {
    const int len1 = static_cast<int>(seq1.length());
    const int len2 = static_cast<int>(seq2.length());

    // Single flat array for the previous row, zero-initialized
    int* prev = new int[len2 + 1]();

    int max_score = 0;

    for (int i = 1; i <= len1; ++i) {
        const char c1 = seq1[i - 1];
        int diag = 0;  // H(i-1, 0) = 0 for all i
        int left = 0;  // H(i,   0) = 0

        for (int j = 1; j <= len2; ++j) {
            const char c2 = seq2[j - 1];

            // Save H(i-1, j) before overwrite
            const int up = prev[j];

            // Diagonal score: match/mismatch
            const int match_score = diag + (c1 == c2 ? MATCH : MISMATCH);

            // Branchless max over all candidates including zero floor
            const int score = std::max({match_score, up + GAP, left + GAP, 0});

            // Track global maximum (branchless via std::max)
            max_score = std::max(max_score, score);

            // Advance diagonal: old H(i-1, j) becomes diag for next column
            diag = up;

            // Store H(i, j) into prev[j] for use by next row
            prev[j] = score;

            // H(i, j) becomes left for next column
            left = score;
        }
    }

    delete[] prev;
    return max_score;
}

struct AlignmentResult {
    int max_score;
    std::string align1;
    std::string align2;
    int max_i, max_j;
};

// Smith-Waterman dynamic programming with Traceback reconstruction
AlignmentResult smith_waterman_traceback(const std::string& seq1, const std::string& seq2) {
    const int len1 = static_cast<int>(seq1.length());
    const int len2 = static_cast<int>(seq2.length());

    // Prevent stack overflow by allocating on the heap
    size_t matrix_size = static_cast<size_t>(len1 + 1) * (len2 + 1);

    // Memory protection limit: ~400 million elements (~1.6 GB RAM)
    // If exceeded, fall back to score-only single-row calculation
    if (matrix_size > 400000000ULL) {
        std::cout << "[Cảnh báo] Kích thước chuỗi lớn. Tự động chuyển sang chế độ Score-only để bảo vệ RAM." << std::endl;
        int score = smith_waterman_sequential(seq1, seq2);
        return {score, "N/A (Chuỗi quá dài)", "N/A (Chuỗi quá dài)", -1, -1};
    }

    int* H = new int[matrix_size]();
    int max_score = 0;
    int max_i = 0, max_j = 0;

    // 1. Fill DP matrix
    for (int i = 1; i <= len1; ++i) {
        const char c1 = seq1[i - 1];
        const size_t row_offset = static_cast<size_t>(i) * (len2 + 1);
        const size_t prev_row_offset = static_cast<size_t>(i - 1) * (len2 + 1);

        for (int j = 1; j <= len2; ++j) {
            const char c2 = seq2[j - 1];

            int score_diag = H[prev_row_offset + (j - 1)] + (c1 == c2 ? MATCH : MISMATCH);
            int score_up   = H[prev_row_offset + j] + GAP;
            int score_left = H[row_offset + (j - 1)] + GAP;

            int score = std::max({score_diag, score_up, score_left, 0});
            H[row_offset + j] = score;

            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // 2. Backtrack to find alignment path
    std::string align1 = "";
    std::string align2 = "";
    int curr_i = max_i;
    int curr_j = max_j;

    while (curr_i > 0 && curr_j > 0) {
        size_t idx = static_cast<size_t>(curr_i) * (len2 + 1) + curr_j;
        int score = H[idx];
        if (score == 0) break; // Stop at local score 0 threshold

        char c1 = seq1[curr_i - 1];
        char c2 = seq2[curr_j - 1];

        int score_diag = H[static_cast<size_t>(curr_i - 1) * (len2 + 1) + (curr_j - 1)];
        int score_up   = H[static_cast<size_t>(curr_i - 1) * (len2 + 1) + curr_j];
        int score_left = H[static_cast<size_t>(curr_i) * (len2 + 1) + (curr_j - 1)];

        if (score == score_diag + (c1 == c2 ? MATCH : MISMATCH)) {
            align1.push_back(c1);
            align2.push_back(c2);
            curr_i--;
            curr_j--;
        } else if (score == score_up + GAP) {
            align1.push_back(c1);
            align2.push_back('-');
            curr_i--;
        } else if (score == score_left + GAP) {
            align1.push_back('-');
            align2.push_back(c2);
            curr_j--;
        } else {
            break;
        }
    }

    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());

    delete[] H;
    return {max_score, align1, align2, max_i, max_j};
}

int main(int argc, char* argv[]) {
    int len1 = 10000;
    int len2 = 10000;
    unsigned int seed = 42;
    bool use_fasta = false;
    std::string file1, file2;

    if (argc >= 3 && is_fasta_file(argv[1]) && is_fasta_file(argv[2])) {
        use_fasta = true;
        file1 = argv[1];
        file2 = argv[2];
    } else {
        if (argc >= 2) len1 = std::atoi(argv[1]);
        if (argc >= 3) len2 = std::atoi(argv[2]);
        if (argc >= 4) seed = static_cast<unsigned int>(std::atoi(argv[3]));
    }

    std::cout << "=== Sequential Smith-Waterman Baseline ===" << std::endl;
    std::string seq1, seq2;

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

    std::cout << "Running Sequential Smith-Waterman..." << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    AlignmentResult result = smith_waterman_traceback(seq1, seq2);
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Max Alignment Score: " << result.max_score << std::endl;
    if (result.max_i != -1) {
        std::cout << "Alignment Ends at Matrix Coordinate: (" << result.max_i << ", " << result.max_j << ")" << std::endl;
        if (result.align1 != "N/A (Chuỗi quá dài)") {
            std::cout << "Reconstructed Alignment Path (First 100 bp or full):" << std::endl;
            size_t print_len = std::min(result.align1.length(), static_cast<size_t>(100));
            std::cout << "  Seq1 (Query): " << result.align1.substr(0, print_len) << (result.align1.length() > 100 ? "..." : "") << std::endl;
            std::cout << "  Seq2 (Db)   : " << result.align2.substr(0, print_len) << (result.align2.length() > 100 ? "..." : "") << std::endl;
            std::cout << "  Alignment Length: " << result.align1.length() << " bp" << std::endl;
        }
    }
    std::cout << "Execution Time: " << elapsed_ms << " ms" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    return 0;
}