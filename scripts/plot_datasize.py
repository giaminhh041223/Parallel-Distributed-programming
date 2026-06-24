import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import sys
import os
import csv

# Set up professional visual aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Segoe UI', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 11,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Path to the sweep results CSV
csv_path = 'sweep_results.csv'

data_sizes = []
total_times = []
comp_times = []

if os.path.exists(csv_path):
    print(f"Loading real sweep benchmark data from {csv_path}")
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_sizes.append(float(row['size']))
            total_times.append(float(row['total_time_s']))
            comp_times.append(float(row['comp_time_s']))
else:
    print("CSV not found — using latest measured cluster data as fallback.")
    data_sizes = [10000, 50000, 100000, 150000, 170000, 180000, 200000]
    total_times = [1.783, 14.2491, 43.6675, 79.1897, 101.680, 125.814, 150.061]
    comp_times = [0.18028, 2.66849, 14.58544, 25.87375, 37.37283, 33.40375, 42.45224]

# Sort data by size for plotting
sorted_idx = np.argsort(data_sizes)
data_sizes = np.array(data_sizes)[sorted_idx]
total_times = np.array(total_times)[sorted_idx]
comp_times = np.array(comp_times)[sorted_idx]

# Plotting
fig, ax = plt.subplots(figsize=(10, 5.5))

# Shade communication overhead first (background layer)
ax.fill_between(data_sizes, comp_times, total_times,
                alpha=0.15, color='#F43F5E', label='Communication & Idle Overhead')

# Plot actual data points and lines (clean professional colors)
ax.plot(data_sizes, total_times, 'o-', color='#0F172A', linewidth=2.0,
        markersize=6, zorder=5, label='Total Time (With Communication)')
ax.plot(data_sizes, comp_times, 's--', color='#0284C7', linewidth=1.8,
        markersize=5, zorder=5, label='Computation Time (Without Communication)')

# Fit quadratic curve (y = a*x^2 + b*x + c)
coeffs = np.polyfit(data_sizes, total_times, 2)
a, b, c = coeffs[0], coeffs[1], coeffs[2]

# Solve a*N^2 + b*N + c = 150 (2.5 mins target)
target_t = 150.0
d = b**2 - 4*a*(c - target_t)
if d >= 0 and a > 0:
    N_target = int((-b + np.sqrt(d)) / (2*a))
else:
    N_target = 300000

# Draw fitted curve (muted color)
x_fit = np.linspace(0, max(data_sizes) * 1.1, 200)
y_fit = np.polyval(coeffs, x_fit)
ax.plot(x_fit, y_fit, ':', color='#8B5CF6', linewidth=1.5, alpha=0.8,
        label=f'Quadratic Fit: {a:.2e}·N² - {abs(b):.2e}·N + {c:.2f}')

# Target zone: 2–3 minutes (120–180 seconds)
ax.axhspan(120, 180, alpha=0.08, color='#10B981', zorder=0)
ax.annotate('Target Zone\n(2–3 min)', xy=(max(x_fit) * 0.95, 150),
            fontsize=9.5, color='#047857', ha='right', va='center',
            fontweight='bold')

# Draw vertical line and plot marker for target N
ax.axvline(x=N_target, color='#10B981', linestyle='--', alpha=0.6, lw=1.2)
ax.plot(N_target, target_t, 'o', color='#10B981', markersize=9, zorder=6,
        markeredgecolor='white', markeredgewidth=1.5,
        label=f'Target N ≈ {N_target:,} bp')

# Labeling & Spines
ax.set_xlabel('Data Size — Length of DNA Sequences (base pairs)', labelpad=8)
ax.set_ylabel('Execution Time (seconds)', labelpad=8)
ax.set_title('Cluster Benchmark vs. Data Size (P = 12 Ranks, Block Size = 128)', fontweight='bold', pad=15)

# Clean grid
ax.grid(True, linestyle='--', alpha=0.5, color='#CBD5E1', zorder=0)

# Remove top and right borders for a modern look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('#94A3B8')
ax.spines['bottom'].set_color('#94A3B8')

ax.legend(loc='upper left', framealpha=0.95, edgecolor='#E2E8F0')
ax.set_xlim(left=0, right=max(data_sizes) * 1.05)
ax.set_ylim(bottom=0, top=max(total_times) * 1.1)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

plt.tight_layout()
plt.savefig('chart_datasize.png', dpi=300)
print("Saved: chart_datasize.png")
print(f"Optimal sequence size N for 2.5 minutes is calculated as: {N_target:,} bp")
