import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import csv

# Set up professional visual aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Segoe UI', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'legend.fontsize': 9.5,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ──────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────
csv_path = sys.argv[1] if len(sys.argv) > 1 else 'benchmark_results.csv'

# Default: use the largest benchmark size
target_size = None
processes   = []
total_time  = []
seq_time    = None

if os.path.exists(csv_path):
    print(f"Loading from {csv_path}")
    data = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            s, p, t = int(row['size']), int(row['num_procs']), float(row['time_ms'])
            data.setdefault(s, {})[p] = t
    target_size = max(data.keys())
    size_data = data[target_size]
    seq_time = size_data.get(0, None)
    for p in sorted(size_data.keys()):
        if p > 0:
            processes.append(p)
            total_time.append(size_data[p])
    print(f"Using size={target_size}, sequential={seq_time:.1f}ms")
else:
    # Default: use the largest benchmark size
    target_size = 600000
    processes  = [1, 2, 4, 12, 24]
    total_time = [4033701.02, 3546819.07, 2202710.00, 1275887.67, 587180.60] # ms
    seq_time   = 4033701.02

# Convert to seconds
total_time_s = [t / 1000.0 for t in total_time]
T1 = (seq_time or total_time[0]) / 1000.0

# ── Derived metrics ──
actual_speedup = [T1 / t for t in total_time_s]
ideal_speedup  = [float(p) for p in processes]
efficiency     = [s / p for s, p in zip(actual_speedup, processes)]

# Estimate Amdahl's f from largest P
S_max = actual_speedup[-1]
P_max = processes[-1]
f_est = (1.0/S_max - 1.0/P_max) / (1.0 - 1.0/P_max)
f_est = max(f_est, 0.001)  # clamp

p_cont = np.linspace(1, P_max * 1.15, 200)
amdahl = 1.0 / (f_est + (1 - f_est) / p_cont)
gustafson = p_cont - f_est * (p_cont - 1)

# ──────────────────────────────────────────────
# 3-subplot figure
# ──────────────────────────────────────────────
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

# Remove top and right spines from all axes
for ax in [ax1, ax2, ax3]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#94A3B8')
    ax.spines['bottom'].set_color('#94A3B8')
    ax.grid(True, ls='--', alpha=0.5, color='#CBD5E1', zorder=0)

# ═══ (a) Execution Time ═══
ax1.plot(processes, total_time_s, 'o-', color='#0F172A', linewidth=2.2, markersize=6,
         label='Measured Wall Time', zorder=5)
ideal_time = [T1 / p for p in processes]
ax1.plot(processes, ideal_time, '--', color='#64748B', linewidth=1.5, alpha=0.8,
         label='Ideal $T_1/P$')

ax1.set_xlabel('Number of Processes (P)', labelpad=8)
ax1.set_ylabel('Execution Time (seconds)', labelpad=8)
ax1.set_title(f'(a) Execution Time (N = {target_size:,})', fontweight='bold', pad=12)
ax1.set_xticks(processes)
ax1.legend(framealpha=0.95, edgecolor='#E2E8F0')
ax1.set_ylim(bottom=0, top=max(total_time_s) * 1.1)

# ═══ (b) Speedup ═══
ax2.plot(processes, actual_speedup, 'o-', color='#0284C7', linewidth=2.2, markersize=6,
         label='Measured Speedup', zorder=5)
ax2.plot(processes, ideal_speedup, '--', color='#64748B', linewidth=1.5, alpha=0.8,
         label='Ideal Linear $S = P$')
ax2.plot(p_cont, amdahl, '-.', color='#8B5CF6', linewidth=1.5, alpha=0.8,
         label=f"Amdahl (f≈{f_est:.4f})")
ax2.plot(p_cont, gustafson, ':', color='#10B981', linewidth=1.5, alpha=0.8,
         label=f"Gustafson (f≈{f_est:.4f})")

# Smart non-overlapping offsets for annotations: (x_offset, y_offset)
offsets = {
    1: (-25, 12),
    2: (12, 10),
    4: (12, -10),
    12: (12, -10),
    24: (-28, -12)
}
for p, s, e in zip(processes, actual_speedup, efficiency):
    ox, oy = offsets[p]
    ax2.annotate(f'{e:.0%}', (p, s),
                 textcoords='offset points', xytext=(ox, oy), fontsize=9,
                 color='#1E293B', fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#94A3B8', lw=0.8))

ax2.set_xlabel('Number of Processes (P)', labelpad=8)
ax2.set_ylabel('Speedup  $S_P = T_1 / T_P$', labelpad=8)
ax2.set_title(f'(b) Speedup Analysis (N = {target_size:,})', fontweight='bold', pad=12)
ax2.set_xticks(processes)
ax2.legend(loc='upper left', framealpha=0.95, edgecolor='#E2E8F0')
ax2.set_ylim(bottom=0, top=max(processes) * 1.1)

# ═══ (c) Parallel Efficiency ═══
# Plot bars evenly spaced by using string categories
processes_str = [str(p) for p in processes]
bar_colors = ['#10B981' if e >= 0.7 else '#F59E0B' if e >= 0.4 else '#EF4444' for e in efficiency]

ax3.bar(processes_str, [e * 100 for e in efficiency], width=0.45,
        color=bar_colors, edgecolor='white', lw=0.6, zorder=3)

# Reference lines
ax3.axhline(100, color='#334155', ls='--', lw=1.2, alpha=0.5, label='100% (Ideal)')
ax3.axhline(70, color='#10B981', ls=':', lw=1.2, alpha=0.6, label='70% (Excellent)')
ax3.axhline(40, color='#F59E0B', ls=':', lw=1.2, alpha=0.6, label='40% (Acceptable)')

for i, e in enumerate(efficiency):
    ax3.text(i, e * 100 + 2, f'{e:.0%}', ha='center', fontsize=9.5,
             fontweight='bold', color='#1E293B')

ax3.set_xlabel('Number of Processes (P)', labelpad=8)
ax3.set_ylabel('Parallel Efficiency $E_P$ (%)', labelpad=8)
ax3.set_title(f'(c) Parallel Efficiency (N = {target_size:,})', fontweight='bold', pad=12)
ax3.set_ylim(0, 115)
ax3.legend(loc='upper right', framealpha=0.95, edgecolor='#E2E8F0')

plt.tight_layout(w_pad=3.0)
plt.savefig('chart_speedup.png', dpi=300)
# plt.show()

# ── Print summary table ──
print(f"\n{'='*55}")
print(f"  SCALABILITY ANALYSIS - N = {target_size:,}")
print(f"{'='*55}")
print(f"{'P':>4s} | {'Time(s)':>9s} | {'Speedup':>8s} | {'Efficiency':>10s}")
print(f"{'-'*4}-+-{'-'*9}-+-{'-'*8}-+-{'-'*10}")
for p, t, s, e in zip(processes, total_time_s, actual_speedup, efficiency):
    print(f"{p:4d} | {t:9.3f} | {s:8.2f} | {e:10.1%}")
print(f"\nAmdahl's sequential fraction f ~= {f_est:.5f}")
print(f"Amdahl's speedup limit S_inf = {1/f_est:.1f}x")
print(f"Gustafson's scaled speedup at P={P_max}: {P_max - f_est*(P_max-1):.1f}x")
print(f"\nSaved: chart_speedup.png")
