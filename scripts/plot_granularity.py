"""
Script 2: Stacked Bar — Load Balance per Process + Granularity Analysis
Purpose: Visualize computation vs. communication/idle per rank.
Reads:   Per-rank timing data (paste from mpi_align output or use sample).
Usage:   python plot_granularity.py
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.titlesize': 15,
    'axes.labelsize': 13,
    'legend.fontsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 250,
    'savefig.bbox': 'tight',
})

# ──────────────────────────────────────────────
# ▶ REPLACE WITH ACTUAL PER-RANK TIMING DATA
#   from mpi_align's [Per-Rank Timing Breakdown]
# ──────────────────────────────────────────────
num_processes = 12

# Computation time per process (seconds) for N=76,700
comp_time = [135.0, 133.5, 134.4, 135.6, 134.1, 134.7,
             135.3, 133.8, 135.0, 135.9, 133.2, 134.4]

# Communication + Idle time per process (seconds) for N=76,700
comm_idle_time = [9.0, 16.5, 15.6, 8.4, 15.9, 15.0,
                  8.7, 16.2, 15.3, 8.1, 16.8, 15.6]

# ──────────────────────────────────────────────

process_ids = list(range(num_processes))
comp_arr  = np.array(comp_time)
comm_arr  = np.array(comm_idle_time)
total_arr = comp_arr + comm_arr

# ── Metrics ──
cv_comp  = np.std(comp_arr) / np.mean(comp_arr) * 100
cv_total = np.std(total_arr) / np.mean(total_arr) * 100
max_imb  = (np.max(total_arr) - np.min(total_arr)) / np.mean(total_arr) * 100
comp_pct = np.mean(comp_arr) / np.mean(total_arr) * 100

# ── Plot ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5),
                                gridspec_kw={'width_ratios': [3, 1]})

# Left: Stacked bar chart
bar_w = 0.55
ax1.bar(process_ids, comp_time, bar_w,
        label='Computation', color='#2471A3', edgecolor='white', lw=0.5)
ax1.bar(process_ids, comm_idle_time, bar_w, bottom=comp_time,
        label='Communication / Idle', color='#E74C3C', edgecolor='white', lw=0.5)

avg_total = np.mean(total_arr)
ax1.axhline(avg_total, color='#2C3E50', ls='--', lw=1.5,
            label=f'Avg Total = {avg_total:.2f}s')

# Annotate % computation on each bar
for i in range(num_processes):
    pct = comp_arr[i] / total_arr[i] * 100
    ax1.text(i, total_arr[i] + 0.05, f'{pct:.0f}%', ha='center',
             fontsize=8, color='#2C3E50', fontweight='bold')

ax1.set_xlabel('Process ID (MPI Rank)')
ax1.set_ylabel('Time (seconds)')
ax1.set_title('Load Distribution Across Processes', fontweight='bold')
ax1.set_xticks(process_ids)
ax1.legend(loc='upper right', framealpha=0.9)
ax1.grid(axis='y', ls='--', alpha=0.3)
ax1.set_ylim(0, max(total_arr) * 1.15)

# Verdict badge
if cv_total < 10:
    verdict, clr = "✅ EXCELLENT", '#27AE60'
elif cv_total < 25:
    verdict, clr = "⚠️  ACCEPTABLE", '#F39C12'
else:
    verdict, clr = "❌ POOR", '#E74C3C'

ax1.text(0.98, 0.97, f'Load Balance: {verdict}',
         transform=ax1.transAxes, fontsize=12, ha='right', va='top',
         color=clr, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=clr, alpha=0.85))

# Right: Summary metrics as a table-like panel
ax2.axis('off')
metrics = [
    ('CV (Computation)', f'{cv_comp:.1f}%'),
    ('CV (Total Time)', f'{cv_total:.1f}%'),
    ('Max Imbalance', f'{max_imb:.1f}%'),
    ('Avg Comp Fraction', f'{comp_pct:.0f}%'),
    ('Avg Comp Time', f'{np.mean(comp_arr):.2f}s'),
    ('Avg Comm/Idle', f'{np.mean(comm_arr):.2f}s'),
    ('Slowest Rank', f'{np.argmax(total_arr)}'),
    ('Fastest Rank', f'{np.argmin(total_arr)}'),
]

ax2.set_title('Summary Metrics', fontweight='bold', fontsize=14)
y_pos = 0.88
for label, value in metrics:
    ax2.text(0.05, y_pos, label + ':', transform=ax2.transAxes,
             fontsize=11, fontweight='bold', color='#2C3E50')
    ax2.text(0.95, y_pos, value, transform=ax2.transAxes,
             fontsize=11, ha='right', family='monospace', color='#34495E')
    y_pos -= 0.10

# Divider
ax2.axhline(y=0.15, xmin=0.05, xmax=0.95, color='#BDC3C7', lw=1,
            transform=ax2.transAxes)

# Interpretation
if cv_total < 10:
    interp = 'Cyclic column distribution\nachieves near-perfect\nload balance across ranks.'
elif cv_total < 25:
    interp = 'Minor imbalance detected.\nConsider tuning Block Size B.'
else:
    interp = 'Significant imbalance!\nReview block distribution\nand network bottlenecks.'

ax2.text(0.50, 0.05, interp, transform=ax2.transAxes,
         fontsize=10, ha='center', va='bottom', style='italic',
         color='#7F8C8D')

plt.tight_layout()
plt.savefig('chart_granularity.png')
plt.show()
print("Saved: chart_granularity.png")
