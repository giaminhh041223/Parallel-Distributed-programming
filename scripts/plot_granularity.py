import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Set up professional visual aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Segoe UI', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'legend.fontsize': 10,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

num_processes = 12

# Raw timing data in milliseconds for N=200,000, B=128 (from physical run)
comp_time_ms = [21624.74, 21418.74, 21557.46, 22329.45, 41325.09, 40243.78, 40576.50, 42452.24, 35130.40, 34893.64, 34804.81, 35879.32]
comm_idle_ms = [128429.15, 128635.26, 128496.55, 127710.96, 108735.64, 109817.20, 109484.77, 107609.20, 114920.48, 115157.32, 115246.22, 114171.79]

# Convert to seconds
comp_time = [t / 1000.0 for t in comp_time_ms]
comm_idle_time = [t / 1000.0 for t in comm_idle_ms]

process_ids = list(range(num_processes))
comp_arr  = np.array(comp_time)
comm_arr  = np.array(comm_idle_time)
total_arr = comp_arr + comm_arr

# ── Metrics ──
cv_comp  = np.std(comp_arr) / np.mean(comp_arr) * 100
cv_total = np.std(total_arr) / np.mean(total_arr) * 100

# Relative idle time deviation according to the formula: (max_idle - min_idle) / max_idle * 100%
max_idle = np.max(comm_arr)
min_idle = np.min(comm_arr)
idle_deviation = (max_idle - min_idle) / max_idle * 100

comp_pct = np.mean(comp_arr) / np.mean(total_arr) * 100

# ── Plot ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                gridspec_kw={'width_ratios': [2.2, 1]})

# Left: Stacked bar chart
bar_w = 0.6
ax1.bar(process_ids, comp_time, bar_w,
        label='Computation', color='#0284C7', edgecolor='white', lw=0.6, zorder=3)
ax1.bar(process_ids, comm_idle_time, bar_w, bottom=comp_time,
        label='Communication / Idle', color='#E11D48', edgecolor='white', lw=0.6, zorder=3)

avg_total = np.mean(total_arr)
ax1.axhline(avg_total, color='#334155', ls='--', lw=1.5, zorder=4,
            label=f'Avg Total = {avg_total:.2f}s')

# Annotate % computation on each bar clearly
for i in range(num_processes):
    pct = comp_arr[i] / total_arr[i] * 100
    ax1.text(i, total_arr[i] + 1.5, f'{pct:.0f}%', ha='center',
             fontsize=8.5, color='#1E293B', fontweight='bold')

ax1.set_xlabel('Process ID (MPI Rank)', labelpad=8)
ax1.set_ylabel('Execution Time (seconds)', labelpad=8)
ax1.set_title('Load Distribution Across Processes (N = 200,000, B = 128)', fontweight='bold', pad=12)
ax1.set_xticks(process_ids)
ax1.legend(loc='upper right', framealpha=0.95, edgecolor='#E2E8F0')
ax1.grid(axis='y', ls='--', alpha=0.5, color='#CBD5E1', zorder=0)
ax1.set_ylim(0, max(total_arr) * 1.15)

# Remove top/right spines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_color('#94A3B8')
ax1.spines['bottom'].set_color('#94A3B8')

# Verdict badge based on teacher's 25% threshold
if idle_deviation <= 25:
    verdict, clr = "ACCEPTABLE", '#10B981'
else:
    verdict, clr = "POOR", '#EF4444'

ax1.text(0.03, 0.95, f'Load Balance: {verdict}',
         transform=ax1.transAxes, fontsize=11, ha='left', va='top',
         color='white', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.4', fc=clr, ec='none', alpha=0.9))

# Right: Summary metrics as a structured table
ax2.axis('off')
table_data = [
    ['Metric', 'Value'],
    ['CV (Computation)', f'{cv_comp:.1f}%'],
    ['Relative Idle Dev', f'{idle_deviation:.1f}%'],
    ['Avg Comp Fraction', f'{comp_pct:.0f}%'],
    ['Avg Comp Time', f'{np.mean(comp_arr):.2f}s'],
    ['Avg Comm/Idle', f'{np.mean(comm_arr):.2f}s'],
    ['Slowest Rank', f'Rank {np.argmax(comp_arr)} (Ryzen)'],
    ['Fastest Rank', f'Rank {np.argmin(comp_arr)} (i7)'],
]

# Create table
tbl = ax2.table(cellText=table_data, loc='center', cellLoc='left')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)
tbl.scale(1.0, 1.8)

# Style cells to look premium
for (row, col), cell in tbl.get_celld().items():
    cell.set_edgecolor('#E2E8F0')
    if row == 0:
        cell.set_text_props(weight='bold', color='#0F172A')
        cell.set_facecolor('#F1F5F9')
    else:
        cell.set_facecolor('white')
        if col == 0:
            cell.set_text_props(color='#475569')
        else:
            cell.set_text_props(color='#0F172A', family='monospace')

ax2.set_title('Performance Metrics', fontweight='bold', pad=15)

# Interpretation text below table
if idle_deviation <= 25:
    interp = f'Idle time deviation is {idle_deviation:.1f}% (<=25%).\nSystem is considered load-balanced.'
else:
    interp = f'Idle deviation is poor ({idle_deviation:.1f}% > 25%).\nWSL2 on Master Node (Core i7)\ncauses computation delay.'

ax2.text(0.50, 0.05, interp, transform=ax2.transAxes,
         fontsize=9.5, ha='center', va='bottom', style='italic',
         color='#64748B')

plt.tight_layout()
plt.savefig('chart_granularity.png', dpi=300)
print(f"Saved: chart_granularity.png with Relative Idle Deviation = {idle_deviation:.2f}%")
