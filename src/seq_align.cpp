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
    int max_score = smith_waterman_sequential(seq1, seq2);
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Max Alignment Score: " << max_score << std::endl;
    std::cout << "Execution Time: " << elapsed_ms << " ms" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    return 0;
}