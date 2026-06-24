# HƯỚNG DẪN SET UP CỤM MPI 3 MÁY (NHÓM BÀI TẬP LỚN)

Tài liệu này ghi lại các bước set up cụm MPI trên 3 laptop Ubuntu của tụi mình để làm bài tập lớn song song. Nhóm mình có 3 người nên cấu hình 1 máy làm Master (máy chính chạy lệnh) và 2 máy làm Slave (máy phụ tính toán).

---

## 1. Chuẩn bị mạng LAN & IP tĩnh

Mấy máy phải cắm chung một cục switch mạng dây (hoặc kết nối chung một mạng Wi-Fi) thì mới thông nhau được. Nhóm mình cần đặt IP tĩnh để lúc chạy không bị đổi IP.

*   **Máy 1 (Master)**: IP là `172.20.10.2` (Máy của bạn A)
*   **Máy 2 (Slave 1)**: IP là `172.20.10.13` (Máy của bạn B - Máy cục bộ)
*   **Máy 3 (Slave 2)**: IP là `172.20.10.3` (Máy của bạn C)

**Test thử**: Mở terminal máy Master gõ ping thử sang 2 máy kia xem có thông không:
```bash
ping 172.20.10.13
ping 172.20.10.3
```
*Lưu ý quan trọng:*
1. Nếu không ping được, nhớ tắt Firewall Windows hoặc chạy lệnh tắt firewall trên Ubuntu: `sudo ufw disable`.
2. Trong môi trường WSL2 kết nối mạng LAN ảo/thực tế, OpenMPI v5 thường gặp lỗi ghép nối giao diện mạng (`Unable to find reachable pairing between local and remote interfaces`). Để khắc phục triệt để lỗi này, **bắt buộc phải vô hiệu hóa IPv6** trên cả 3 máy bằng các lệnh sau:
   ```bash
   sudo sysctl -w net.ipv6.conf.all.disable_ipv6=1
   sudo sysctl -w net.ipv6.conf.default.disable_ipv6=1
   ```

---

## 2. Đồng bộ Username (Cực kỳ quan trọng)

MPI nó bắt buộc username trên cả 3 máy phải giống hệt nhau thì SSH mới tự động đăng nhập được.
Nhóm mình thống nhất tạo tài khoản chung tên là `hpcuser` trên cả 3 máy:
```bash
# Tạo user mới
sudo adduser hpcuser
# Cấp quyền admin
sudo usermod -aG sudo hpcuser
```
*Sau khi tạo xong thì nhớ Log out tài khoản cũ rồi Log in vào tài khoản `hpcuser` mới này trên cả 3 máy để làm tiếp.*

---

## 3. Cài đặt các gói phần mềm cần thiết

Chạy lệnh này trên **cả 3 máy** (mở terminal lên gõ):
```bash
sudo apt update
sudo apt install -y build-essential openmpi-bin openmpi-common libopenmpi-dev openssh-server openssh-client rsync
```

---

## 4. Cấu hình SSH Key để kết nối không cần mật khẩu

Để máy Master tự động nhảy sang các máy Slave chạy file mà không bị hỏi mật khẩu:

1.  **Trên máy Master (172.20.10.2)**, chạy lệnh sinh key:
    ```bash
    # Cứ nhấn Enter liên tiếp cho qua, không nhập passphrase
    ssh-keygen -t rsa
    ```
2.  **Copy key từ Master sang 2 máy Slave và chính nó**:
    ```bash
    ssh-copy-id hpcuser@172.20.10.2
    ssh-copy-id hpcuser@172.20.10.13
    ssh-copy-id hpcuser@172.20.10.3
    ```
3.  **Test thử xem được chưa**: Đứng ở Master gõ `ssh hpcuser@172.20.10.13`. Nếu nó đăng nhập thẳng vào luôn không đòi mật khẩu là ok. Nhớ gõ `exit` để thoát ra nhé.

---

## 5. Đồng bộ thư mục làm việc

File chạy `mpi_align` phải nằm ở cùng một đường dẫn trên cả 3 máy thì MPI mới tìm thấy để chạy.
Tụi mình thống nhất để ở đường dẫn: `/home/hpcuser/smith_waterman`

1.  Tạo thư mục này trên cả 3 máy:
    ```bash
    mkdir -p /home/hpcuser/smith_waterman
    ```
2.  Biên dịch code trên máy Master:
    ```bash
    cd /home/hpcuser/smith_waterman
    mpicxx -O3 -o mpi_align mpi_align.cpp
    ```
3.  Đồng bộ file chạy sang 2 máy Slave (mỗi lần sửa code xong biên dịch lại thì chạy lệnh này để copy sang các máy kia):
    ```bash
    rsync -avz /home/hpcuser/smith_waterman/ hpcuser@172.20.10.13:/home/hpcuser/smith_waterman/
    rsync -avz /home/hpcuser/smith_waterman/ hpcuser@172.20.10.3:/home/hpcuser/smith_waterman/
    ```

---

## 6. Chạy thử chương trình

1.  **Tạo file hostfile**: Trên máy Master, tạo một file text tên là `hosts.cfg` nằm trong `/home/hpcuser/smith_waterman/`. Nội dung điền IP và số nhân CPU (slots) của mỗi máy:
    ```text
    172.20.10.2 slots=4
    172.20.10.13 slots=4
    172.20.10.3 slots=4
    ```
2.  **Gõ lệnh chạy song song**:
    ```bash
    cd /home/hpcuser/smith_waterman
    # Chạy thử với 12 tiến trình (mỗi máy 4 cores)
    mpirun -np 12 --hostfile hosts.cfg ./mpi_align 76700 76700 42 128
    ```
    Hoặc chạy kịch bản tự động test bằng file bash:
    ```bash
    ./run_test.sh hosts.cfg
    ```
3.  **Chạy thử nghiệm với dữ liệu sinh học thực tế (SARS-CoV-2 vs MERS-CoV)**:
    - **Bước 1**: Đảm bảo thư mục `/home/hpcuser/smith_waterman/data/` trên máy Master chứa các file `sars_cov_2.fasta` và `mers_cov.fasta`. (Lưu ý: Chỉ cần máy Master có các file này là đủ, các máy Slave sẽ nhận dữ liệu qua mạng LAN bằng lệnh phát sóng MPI_Bcast).
    - **Bước 2**: Chạy demo song song kết hợp so sánh tuần tự bằng file script:
      ```bash
      chmod +x run_demo_fasta.sh
      ./run_demo_fasta.sh hosts.cfg
      ```
      Kết quả điểm so khớp (**Max Score**) in ra phải trùng khớp hoàn toàn là **30789** giữa cả hai phiên bản.
