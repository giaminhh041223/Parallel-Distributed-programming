# Parallel Smith-Waterman Sequence Alignment (OpenMPI)

## 📌 Tổng quan dự án (Project Overview)
Dự án song song hóa thuật toán **Smith-Waterman (Local Sequence Alignment)** phục vụ cho việc so khớp chuỗi DNA sinh học cục bộ. Thuật toán được tối ưu hóa ở cấp độ HPC (High-Performance Computing) và song song hóa trên hệ thống **bộ nhớ phân tán (Distributed Memory)** sử dụng thư viện **OpenMPI**.

Dự án hỗ trợ chạy kiểm thử cả trên bộ dữ liệu sinh ngẫu nhiên có kiểm soát lẫn **dữ liệu sinh học thực tế** (Hệ gen virus SARS-CoV-2 vs MERS-CoV) được tải trực tiếp từ ngân hàng gen quốc gia Hoa Kỳ **NCBI**.

---

## 🌟 Các tính năng nổi bật (Key Features)
1.  **Tối ưu hóa tuần tự (Sequential Optimization)**:
    *   Sử dụng bộ đệm cuộn 1 dòng (**Single-row rolling buffer**) giúp giảm độ phức tạp không gian xuống $O(N)$ thay vì ma trận $O(M \times N)$ khổng lồ.
    *   Áp dụng tính toán không rẽ nhánh (**Branchless std::max**) giúp CPU phát ra lệnh dịch chuyển có điều kiện (`cmov`) thay vì rẽ nhánh gây chậm trễ.
2.  **Tối ưu hóa song song MPI (OpenMPI Wavefront)**:
    *   Chia khối 2D (**2D Block-based Decomposition**) kết hợp tính toán chéo Wavefront Pipeline.
    *   Phân phối cột khối xoay vòng (**Cyclic Column Block assignment**) để tự động tối ưu hóa cân bằng tải (Load Balancing) giữa các máy.
    *   Nạp trước bộ nhớ cache (**Memory Prefetching**) tương thích đa nền tảng để tránh trễ truy xuất RAM.
    *   Giao tiếp không chặn hợp nhất (**Consolidated Non-blocking communication**) tránh deadlock mạng LAN.
3.  **Tích hợp dữ liệu y sinh NCBI**:
    *   Đọc và xử lý định dạng file sinh học **FASTA** trực tiếp.
    *   Hỗ trợ cơ chế Broadcast hiệu quả: Chỉ cần máy Master đọc file FASTA từ đĩa cứng, dữ liệu sẽ được phát sóng tự động tới RAM của các máy Slave qua mạng LAN mà không cần copy file vật lý sang các máy con.

---

## 📁 Cấu trúc thư mục (Directory Structure)
*   [src/seq_align.cpp](src/seq_align.cpp): Phiên bản thuật toán tuần tự đã tối ưu hóa.
*   [src/mpi_align.cpp](src/mpi_align.cpp): Phiên bản song song hóa MPI.
*   [data/](data/): Thư mục chứa các file hệ gen của SARS-CoV-2, SARS-CoV-1, MERS-CoV dạng FASTA.
*   [scripts/](scripts/): Các script vẽ biểu đồ hiệu năng (Speedup, Granularity) bằng Python.
*   [run_demo_fasta_windows.bat](run_demo_fasta_windows.bat): Script chạy demo tự động so sánh tuần tự và song song trên Windows.
*   [run_demo_fasta.sh](run_demo_fasta.sh): Script chạy demo tự động trên cụm Linux Master/Slaves.
*   [run_test.sh](run_test.sh): Script kiểm thử độ chính xác và benchmark quét các kích thước ma trận.

---

## 🚀 Hướng dẫn chạy chương trình (How to Run)

### 💻 Trên Windows (Chạy cục bộ giả lập)
Yêu cầu hệ thống đã cài đặt **Visual Studio (C++)** và **MS-MPI**.
1.  Mở thư mục dự án.
2.  Nhấp đúp chuột vào file [run_demo_fasta_windows.bat](run_demo_fasta_windows.bat) để khởi chạy chạy demo so khớp SARS-CoV-2 vs MERS-CoV.

### 🐧 Trên cụm máy Ubuntu vật lý (3 nodes - 12 cores)
Yêu cầu cả 3 máy đã cài đặt **OpenMPI** và thông kết nối **SSH không mật khẩu** thông qua tài khoản chung `hpcuser`.

1.  Tạo file `hosts.cfg` trên máy Master chỉ định IP và số nhân CPU khả dụng:
    ```text
    192.168.1.100 slots=4
    192.168.1.101 slots=4
    192.168.1.102 slots=4
    ```
2.  Chạy script demo so khớp gen sinh học thực tế:
    ```bash
    chmod +x run_demo_fasta.sh
    ./run_demo_fasta.sh hosts.cfg
    ```
3.  Kết quả hiển thị điểm số tối ưu và bảng đo đạc thời gian phân rã của từng máy tính.

---

## 📈 Kết quả thử nghiệm tiêu biểu (Sample Benchmarks)
So khớp dữ liệu thật: **SARS-CoV-2 (29,903 bp) vs MERS-CoV (30,119 bp)**
*   **Max Alignment Score**: **`30789`** (Khớp hoàn hảo 100% giữa hai phiên bản)
*   **Sequential Baseline**: `20.26` giây.
*   **OpenMPI (3 Ranks)**: `14.04` giây.
*   **Hệ số tăng tốc (Speedup)**: **`1.44x`** (đo trên môi trường Windows loopback).
