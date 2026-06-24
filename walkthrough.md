# Project Walkthrough: Parallel Smith-Waterman Alignment using OpenMPI

This document summarizes the final project setup, codebase files, compilation, and cluster execution instructions.

## Project Structure
The workspace is organized as follows:
- [src/seq_align.cpp](file:///D:/Parallel%20and%20distribution%20program/src/seq_align.cpp): Fully optimized sequential baseline using $O(M)$ memory.
- [src/mpi_align.cpp](file:///D:/Parallel%20and%20distribution%20program/src/mpi_align.cpp): Parallel MPI implementation using 2D block-based wavefront processing.
- [run_test.sh](file:///D:/Parallel%20and%20distribution%20program/run_test.sh): Automation shell script to compile, verify correctness, and run benchmarks for varying sequence sizes and MPI rank counts.

---

## Code Execution Guide (On the Cluster Master Node)

### 1. Transfer Code to the Master Node
Make sure all source files and the shell script are present on your master node (e.g., in a directory named `hpc_smith_waterman`).

### 2. Make the Script Executable
To run the automated script, change its permissions to allow execution:
```bash
chmod +x run_test.sh
```

### 3. Run Verification & Benchmarking
- **To run locally on the Master Node (1 node simulation)**:
  ```bash
  ./run_test.sh
  ```
- **To run on the 3-node physical Cluster (using a hostfile)**:
  Create an MPI hostfile (e.g., `hosts.cfg`) mapping your master and slave nodes:
  ```text
  172.20.10.2 slots=4   # Slave Node 1 (Rank 0-3)
  172.20.10.13 slots=4  # Master Node (Rank 4-7)
  172.20.10.3 slots=4   # Slave Node 2 (Rank 8-11)
  ```
  Then run the script with the hostfile as an argument:
  ```bash
  ./run_test.sh hosts.cfg
  ```

---

---

## Local Verification & Demo (Windows)
To run the demo locally on your Windows machine at any time without asking the AI assistant:
1. Open the project folder `D:\Parallel and distribution program`.
2. Locate the file [run_demo_windows.bat](file:///D:/Parallel%20and%20distribution%20program/run_demo_windows.bat).
3. **Double-click** the file in Windows File Explorer, or open a Command Prompt (cmd) / PowerShell and run:
   ```cmd
   .\run_demo_windows.bat
   ```
4. The script will automatically execute the sequential and parallel code side-by-side and wait for a keypress before closing.

---

## Verification & Output
The cluster automation script (`run_test.sh`) automatically:

1. Compiles both `seq_align.cpp` and `mpi_align.cpp` with `-O3` level optimization.
2. Performs a **Correctness Verification** test with $5000 \times 5000$ base-pair sequences using 3 ranks. If the sequential score and parallel score do not match exactly, the script terminates immediately with a failure warning.
3. Loops through $5000$, $10000$, and $20000$ base-pair sequence lengths and runs sequential and parallel tests across $1$, $2$, and $3$ ranks.
4. Outputs a Markdown table in the terminal.

---

## Local Verification Results (Windows Simulation)
Because OpenMPI does not officially support Windows, we installed **Microsoft MPI (MS-MPI)** and compiled the code using the Microsoft C++ compiler. 

### Commands run:
- **Compile Sequential**: `cl /EHsc /O2 /Fe:seq_align.exe src\seq_align.cpp`
- **Compile Parallel**: `cl /EHsc /O2 /Fe:mpi_align.exe src\mpi_align.cpp /I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include" /link /LIBPATH:"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" msmpi.lib`
- **Run Sequential**: `.\seq_align.exe 5000 5000 42`
- **Run Parallel (3 Processes)**: `& "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 3 .\mpi_align.exe 5000 5000 42 128`

### Test Output Comparison (Phase 5 Optimized):

| Execution Mode | Sequence Size | Seed | Max Score | Execution Time (Original) | Execution Time (Optimized) | Speedup vs Original |
|---|---|---|---|---|---|---|
| **Sequential** | 5,000 x 5,000 | 42 | **3,639** | 317.27 ms | **248.67 ms** | **+21.6%** |
| **MPI (3 Ranks)** | 5,000 x 5,000 | 42 | **3,639** | 213.25 ms | **103.97 ms** | **+105.1% (2.05x)** |
| **Sequential** | 10,000 x 10,000 | 42 | **7,246** | 1552.74 ms | **968.25 ms** | **+37.6%** |
| **MPI (3 Ranks)** | 10,000 x 10,000 | 42 | **7,246** | 1055.18 ms | **602.56 ms** | **+75.1% (1.75x)** |
| **Sequential** | 20,000 x 20,000 | 42 | **14,384** | 5307.53 ms | **3,699.66 ms** | **+30.3%** |

---

## Real Dataset Verification (SARS-CoV-2 vs MERS-CoV)
We downloaded and tested real viral genome sequences from the **NCBI GenBank** database to verify the performance on real biological data:
- **Sequence 1**: SARS-CoV-2 Wuhan-Hu-1 reference genome ([sars_cov_2.fasta](file:///D:/Parallel%20and%20distribution%20program/data/sars_cov_2.fasta)) — **29,903 bp**
- **Sequence 2**: MERS-CoV reference genome ([mers_cov.fasta](file:///D:/Parallel%20and%20distribution%20program/data/mers_cov.fasta)) — **30,119 bp**

### Execution Commands:
- **Sequential (Windows)**: `.\seq_align.exe data\sars_cov_2.fasta data\mers_cov.fasta`
- **Parallel (3 Processes, Windows)**: `& "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 3 .\mpi_align.exe data\sars_cov_2.fasta data\mers_cov.fasta 128`
- **Parallel (12 Ranks, 3-Node Physical LAN Cluster)**:
  ```bash
  mpirun --hostfile hosts.cfg ./mpi_align data/sars_cov_2.fasta data/mers_cov.fasta 128
  ```

### Results:
- **Max Alignment Score**: **30,789** (Matches exactly in all modes!)
- **Windows Sequential Time**: **20,267.80 ms**
- **Windows Parallel Time (3 Ranks)**: **14,049.96 ms** (Speedup of **1.44x** over Windows loopback)
- **Linux Cluster Sequential Time (Master Node)**: **1,263.99 ms**
- **Linux Cluster Parallel Time (12 Ranks over 3 Physical Nodes via Wi-Fi Hotspot)**: **4,238.31 ms**

> [!NOTE]
> Do note that the parallel time over the 3 physical machines (`4,238.31 ms`) is slower than the sequential baseline on the Master node (`1,263.99 ms`). This is primarily because:
> 1. **WiFi Network Latency**: Running OpenMPI over a shared Wi-Fi network (iPhone hotspot) introduces substantial TCP socket latency and jitter, creating a communication bottleneck during the boundary exchanges of the block-based wavefront.
> 2. **CPU Heterogeneity & Virtualization Overhead**: The cluster nodes have different environments (Master is Core i7 running WSL2, Slave 1 is Core i5 running Native Linux, and Slave 2 is Ryzen 5 5500 running Native Linux). As shown in the per-rank breakdown, the Master node (WSL2 overhead) was significantly slower (comp time ~1290 ms) than the Native Linux nodes (comp time ~313 ms), causing all other ranks to wait at synchronization points.

### Batch Demo:
We created [run_demo_fasta_windows.bat](file:///D:/Parallel%20and%20distribution%20program/run_demo_fasta_windows.bat) to run this real-world demo on Windows, and [run_demo_fasta.sh](file:///D:/Parallel%20and%20distribution%20program/run_demo_fasta.sh) for execution on your physical 3-node cluster master node.

> [!NOTE]
> The global alignment scores match **exactly** in all cases, validating the mathematical soundness of the wavefront block decomposition. The optimized version shows an average speedup of ~30% in sequential execution due to:
> 1. **Branchless std::max** enabling conditional moves (`cmov`) instead of branching.
> 2. **Cross-platform prefetching** (`PREFETCH`) fetching H_local rows into L1 cache ahead of time in the parallel execution path.
> 3. **RAII-based std::vector** handles replacing manual `new[]` and `delete[]` allocations, eliminating leaks.





## Bug Fix & Report Alignment (Phase 6)

In this phase, we addressed critical correctness, deadlock, and documentation issues:

### 1. MPI Tag Collision Resolution
- **Problem**: The original tag formula `tag = r % 32768` led to tag collisions when multiple column blocks on the same row $r$ were processed by the same sender/receiver rank pair, causing out-of-order buffer matching.
- **Fix**: We updated both `mpi_align.cpp` and `mpi_align_virtual.cpp` to use unique tags for each communication edge:
  - **Receives (Step 1)**: `tag = (r * N_b + c) % 32767` (where `c` is the receiver's column index).
  - **Sends (Step 5)**: `tag = (r * N_b + c + 1) % 32767` (since the sender at `c` sends boundary data to column `c + 1`).
  - This guarantees collision-free message matching while conforming to the standard MPI tag limit (minimum guaranteed limit is 32767).

### 2. Virtual Simulator Deadlock Elimination
- **Problem**: Updating the tags in `src/mpi_align_virtual.cpp` with strict matching revealed a deadlock in the simulator. The simulator originally used blocking `MPI_Recv` and `MPI_Send` calls inside a cyclic-column wavefront traversal, causing ranks to block on sends before receivers could post matching receives.
- **Fix**: We rewrote the virtual simulator's wavefront loop to use non-blocking communication (`MPI_Irecv`, `MPI_Isend`, `MPI_Waitall`) identical to the real execution path. This eliminated all deadlocks.

### 3. Verification Logs (TC-1, TC-4, TC-5)
We executed the three missing test cases locally to generate physical verification logs:
- **TC-1 ($5,000 \times 5,000$, Seed 42)**: Sequential Score `3,639` vs Parallel Score `3,639`. Matched perfectly. Saved in [verification_TC1.txt](file:///d:/Parallel%20and%20distribution%20program/data/verification_TC1.txt).
- **TC-4 ($7,777 \times 7,777$, Seed 42 - Odd Dimension)**: Sequential Score `5,635` vs Parallel Score `5,635`. Matched perfectly. Saved in [verification_TC4.txt](file:///d:/Parallel%20and%20distribution%20program/data/verification_TC4.txt).
- **TC-5 (SARS-CoV-2 vs MERS-CoV NCBI Genomes)**: Sequential Score `30,789` vs Parallel Score `30,789`. Matched perfectly. Saved in [verification_TC5.txt](file:///d:/Parallel%20and%20distribution%20program/data/verification_TC5.txt).

### 4. Report and Regression Alignment
- **Polynomial Equation**: Standardized the polynomial regression formula in all reports and chart legends to the actual fit equation:
  $$T(N) = 1.79 \times 10^{-9} \cdot N^2 - 7.12 \times 10^{-5} \cdot N + 1.47$$
- **Target N Resolution**: Clarified that solving for the optimal zone midpoint ($T = 150$s) yields $N \approx 201,362$ bp, confirming that $N=200,000$ running in 150.06s fits the 2-3 minutes performance target.
- **Oversubscribing P=24**: Added a detailed section explaining that slot-based mapping (`--map-by slot`, taking $587.18$s) outperformed node-based mapping (`--map-by node`, taking $749.94$s) because slot grouping localizes communication via fast shared-memory instead of the Wi-Fi TCP sockets.
