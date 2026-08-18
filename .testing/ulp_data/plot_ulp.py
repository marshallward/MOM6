#!/usr/bin/env python3
"""Plot ULP error vs x for exp_repro and intrinsic exp()."""

import matplotlib.pyplot as plt
import csv

def load_csv(filename):
    x = []
    ulp = []
    with open(filename) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            x.append(float(row[0]))
            ulp.append(float(row[2]))
    return x, ulp

# Load data
repro_scalar_x, repro_scalar_ulp = load_csv('exp_repro_scalar_ulp.csv')
repro_vector_x, repro_vector_ulp = load_csv('exp_repro_vector_ulp.csv')
exp_scalar_x, exp_scalar_ulp = load_csv('exp_scalar_ulp.csv')
exp_vector_x, exp_vector_ulp = load_csv('exp_vector_ulp.csv')

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)

# exp_repro scalar
ax = axes[0, 0]
ax.scatter(repro_scalar_x, repro_scalar_ulp, s=1, alpha=0.5)
ax.set_ylabel('ULP error')
ax.set_title('exp_repro (scalar)')
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.5)
ax.axhline(y=1.0, color='gray', linestyle='-', linewidth=0.5)

# exp_repro vector
ax = axes[0, 1]
ax.scatter(repro_vector_x, repro_vector_ulp, s=1, alpha=0.5)
ax.set_title('exp_repro (vector)')
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.5)
ax.axhline(y=1.0, color='gray', linestyle='-', linewidth=0.5)

# exp() scalar
ax = axes[1, 0]
ax.scatter(exp_scalar_x, exp_scalar_ulp, s=1, alpha=0.5, color='C1')
ax.set_xlabel('x')
ax.set_ylabel('ULP error')
ax.set_title('exp() intrinsic (scalar)')
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.5)
ax.axhline(y=1.0, color='gray', linestyle='-', linewidth=0.5)

# exp() vector
ax = axes[1, 1]
ax.scatter(exp_vector_x, exp_vector_ulp, s=1, alpha=0.5, color='C1')
ax.set_xlabel('x')
ax.set_title('exp() intrinsic (vector)')
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.5)
ax.axhline(y=1.0, color='gray', linestyle='-', linewidth=0.5)

for ax in axes.flat:
    ax.set_ylim(0, 1.5)

plt.suptitle('ULP error vs x [-10, 10] (GCC)')
plt.tight_layout()
plt.savefig('ulp_comparison.png', dpi=150)
print('Saved ulp_comparison.png')
