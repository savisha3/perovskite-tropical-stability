"""Re-run the stress budget for all 20 sites and compare with Table 1 of the review.

Usage (from the repository root):
  python tests/reproduce_table1.py
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import stressbudget as sb  # noqa: E402

COLS = {
    'damp_heat_open_Ea0.6': 'open damp heat 1000 h, Ea 0.6',
    'damp_heat_open_Ea0.3': 'open damp heat 1000 h, Ea 0.3',
    'damp_heat_open_Ea0.95': 'open damp heat 1000 h, Ea 0.95',
    'light_heat_open_Ea0.61': 'open light+heat 85C 1000 h, Ea 0.61',
    'light_heat_open_Ea0.43': 'open light+heat 85C 1000 h, Ea 0.43',
    'dark_heat_open': 'open dark heat 85C 1000 h',
    'thermal_cycling_open': 'open thermal cycling 200 cycles',
    'thermal_cycling_roof': 'roof thermal cycling 200 cycles',
}

cities = pd.read_csv(os.path.join(ROOT, 'data', 'cities.csv'))
ref = pd.read_csv(os.path.join(ROOT, 'tests', 'table1_reference.csv'), index_col='city')
worst = 0.0
bad = 0
for _, c in cities.iterrows():
    s = sb.summary(os.path.join(ROOT, c.file), c.tilt)
    for rc, sc in COLS.items():
        r, v = ref.loc[c.city, rc], s[sc]
        rel = abs(v - r) / r
        worst = max(worst, rel)
        if abs(v - r) > 0.0051 + 1e-9 and rel > 0.001:
            bad += 1
            print(f'MISMATCH {c.city} {rc}: reference {r}, now {v:.3f}')
    print(f'{c.city:17s} damp heat {s[COLS["damp_heat_open_Ea0.6"]]:6.2f} yr   light+heat {s[COLS["light_heat_open_Ea0.61"]]:5.2f} yr')
print(f'Largest relative difference from Table 1 (reference values rounded to 0.01 yr): {worst:.2e}')
print('ALL MATCH' if bad == 0 else f'{bad} MISMATCHES')
