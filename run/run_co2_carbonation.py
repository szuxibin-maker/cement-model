#!/usr/bin/env python3
"""Standalone GEMS CO₂ carbonation simulation script.

Simulates the effect of progressively increasing CO₂ amounts on a fully
hydrated cement paste and generates results similar to Figure 14:

  (a) Stacked area plot of phase volume fractions + pH overlay curve
  (b) Semi-log plot of ionic concentrations vs CO₂ amount

Fixed cement recipe (CEM III/A blend)
--------------------------------------
  Clinkers : C3S, C2S, C3A, C4AF  (from Bogue's composition)
  SCM      : GGBFS
  Gypsum   : Gp  (CaSO₄·2H₂O)
  Water    : H2O@
  Oxygen   : O2  (trace — reduces numerical stiffness)

Usage (run from project root)
------------------------------
  python run/run_co2_carbonation.py

Output files (written to the same directory as this script)
------------------------------------------------------------
  co2_results.csv
  Figure_14a_phase_evolution.png
  Figure_14b_ionic_concentration.png
"""

import os
import sys

# ---------------------------------------------------------------------------
# Make sure the project root is on sys.path so we can import run.* / util.*
# ---------------------------------------------------------------------------
_BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _BASE_DIR not in sys.path:
    sys.path.insert(0, _BASE_DIR)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')                     # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from run.GEMSCalc import GEMS
from util.final_hydration import SCM      # elemental formulas for SCMs

# ===========================================================================
# Configuration
# ===========================================================================

# Cement recipe  (g per 100 g cement blend)
CLINK_PHASES = {
    "C3S":  50.0,   # Alite
    "C2S":  15.0,   # Belite
    "C3A":   8.0,   # Aluminate
    "C4AF":  8.0,   # Ferrite
}
GGBFS_CONTENT  = 50.0   # g / 100 g blend
GYPSUM_CONTENT =  5.0   # g / 100 g blend  (Gp = CaSO₄·2H₂O)
WATER_CONTENT  = 50.0   # g / 100 g blend  (w/b ≈ 0.50)
TEMPERATURE    = 25.0   # °C

# CO₂ scan
CO2_MIN   =  0.0   # g / 100 g blend
CO2_MAX   = 60.0
CO2_STEPS = 13     # 0, 5, 10, …, 60

# GEMS database
_GEMS_FILE = os.path.join(_BASE_DIR, 'gems_files', 'MySystem2-dat.lst')

# Output paths
_OUT_DIR      = os.path.dirname(os.path.abspath(__file__))
CSV_OUTPUT    = os.path.join(_OUT_DIR, 'co2_results.csv')
FIG_PHASE_OUT = os.path.join(_OUT_DIR, 'Figure_14a_phase_evolution.png')
FIG_IONIC_OUT = os.path.join(_OUT_DIR, 'Figure_14b_ionic_concentration.png')

# Ionic species to monitor: GEMS species name → display label
IONIC_SPECIES = {
    'Ca+2':    r'Ca$^{2+}$',
    'K+':      r'K$^+$',
    'Na+':     r'Na$^+$',
    'Al+3':    r'Al$^{3+}$',
    'H4SiO4@': r'Si(OH)$_4$',
    'OH-':     r'OH$^-$',
    'Mg+2':    r'Mg$^{2+}$',
    'SO4-2':   r'SO$_4^{2-}$',
}

# ===========================================================================
# Phase suppression  (consistent with util/final_hydration.py)
# ===========================================================================

def _suppress_phases(gemsk):
    """Suppress phases irrelevant for OPC + GGBFS carbonation."""
    _to_suppress = [
        # Rare / high-temperature AFt variants
        "thaumasite",
        # Chloride-bearing phases
        "Friedels",
        # Zeolites (synthetic / not natural cement products)
        "ZeoliteP", "Natrolite", "Chabazite", "ZeoliteX", "ZeoliteY",
        # High-alumina clinker phases
        "CA", "CA2", "Mayenite",
        # Free lime and alkali salts
        "lime", "arcanite", "thenardite", "Na-oxide", "K-oxide",
        # CaCO₃ metastable polymorph
        "Aragonite",
        # Iron minerals (Fe₂O₃ too low in OPC to form these)
        "Iron", "Fe-carbonate", "Siderite", "Hematite", "Magnetite",
        "Ferrihydrite-am", "Ferrihydrite-mc", "Goethite", "Pyrite",
        "Troilite", "Melanterite",
        # Manganese minerals
        "Rhodochrosite", "Rhodochrosite-sy", "Hausmannite",
        "Pyrolusite", "Manganite", "Pyrochroite",
        # Misc.
        "Kaolinite", "Graphite",
        "Dolomite-dis", "Dolomite-ord",
        "syngenite",
        "Ti(alpha)", "TiO2(am_hyd)",
    ]
    for name in _to_suppress:
        try:
            gemsk.supress_phase(name)
        except Exception:
            pass  # Phase may not exist in this database; skip silently

# ===========================================================================
# Core simulation
# ===========================================================================

def _setup_gems(co2_g):
    """Return a configured GEMS instance for *co2_g* g of CO₂ / 100 g blend."""
    gemsk = GEMS(_GEMS_FILE)
    _suppress_phases(gemsk)

    # Clinker phases
    for phase, amount in CLINK_PHASES.items():
        gemsk.add_species_amt(phase, amount * 1e-3, units="kg")

    # Gypsum  (Gp = CaSO₄·2H₂O)
    gemsk.add_species_amt("Gp", GYPSUM_CONTENT * 1e-3, units="kg")

    # GGBFS — added as elemental formula (from util/final_hydration.py SCM dict)
    ggbfs_formula, ggbfs_units = SCM['GGBFS']
    gemsk.add_amt_from_formula(ggbfs_formula, GGBFS_CONTENT * 1e-3, units=ggbfs_units)

    # Water
    gemsk.add_species_amt("H2O@", WATER_CONTENT * 1e-3, units="kg")

    # Trace O₂ to reduce numerical stiffness (redox)
    gemsk.add_species_amt("O2", 1e-6)

    # CO₂
    if co2_g > 0.0:
        gemsk.add_species_amt("CO2", co2_g * 1e-3, units="kg")

    gemsk.T = TEMPERATURE + 273.15
    return gemsk


def run_co2_scan():
    """
    Run GEMS equilibration for each CO₂ level.

    Returns
    -------
    co2_values : np.ndarray
    results    : list of dict (or None when equilibration failed)
    """
    co2_values = np.linspace(CO2_MIN, CO2_MAX, CO2_STEPS)
    results = []

    print(f"CO₂ carbonation scan  —  {CO2_STEPS} steps  "
          f"[{CO2_MIN:.0f} … {CO2_MAX:.0f} g / 100 g blend]")
    print("-" * 65)

    for idx, co2 in enumerate(co2_values):
        print(f"  Step {idx + 1:>2}/{CO2_STEPS}  CO₂ = {co2:5.1f} g … ",
              end="", flush=True)
        try:
            gemsk = _setup_gems(co2)
            status = gemsk.equilibrate()

            if status.startswith('F') or status.startswith('B'):
                print(f"FAILED  [{status}]")
                results.append(None)
                continue

            # Phase volume fractions (gas phase excluded from denominator)
            phase_vols = gemsk.phase_volumes
            total_vol  = sum(v for k, v in phase_vols.items() if k != 'gas_gen')
            phase_vfrac = {
                k: (v / total_vol if total_vol > 0 else 0.0)
                for k, v in phase_vols.items()
                if k != 'gas_gen'
            }

            results.append({
                'co2':            co2,
                'status':         status,
                'pH':             gemsk.pH,
                'pE':             gemsk.pE,
                'ionic_strength': gemsk.ionic_strength,
                'phase_volumes':  phase_vols,
                'phase_masses':   gemsk.phase_masses,
                'phase_vfrac':    phase_vfrac,
                'aq_composition': gemsk.aq_composition,
            })
            print(f"OK  (pH = {results[-1]['pH']:.2f})")

        except Exception as exc:
            print(f"ERROR — {exc}")
            results.append(None)

    return co2_values, results

# ===========================================================================
# Data export
# ===========================================================================

def export_csv(co2_values, results):
    """Save all simulation data to *CSV_OUTPUT* and return the DataFrame."""
    rows = []
    for co2, res in zip(co2_values, results):
        if res is None:
            continue
        row = {
            'co2_g_per_100g': co2,
            'pH':             res['pH'],
            'pE':             res['pE'],
            'ionic_strength': res['ionic_strength'],
        }
        for phase, vfrac in res['phase_vfrac'].items():
            row[f'vfrac_{phase}'] = vfrac
        for species, conc in res['aq_composition'].items():
            row[f'aq_{species}'] = conc
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(CSV_OUTPUT, index=False)
    print(f"\nCSV saved  →  {CSV_OUTPUT}")
    return df

# ===========================================================================
# Plotting
# ===========================================================================

def plot_phase_evolution(co2_values, results, save_path=None):
    """
    Figure 14(a) — stacked area of phase volume fractions with pH overlay.

    Parameters
    ----------
    co2_values : array-like
    results    : list of dict (None entries are skipped)
    save_path  : str or None  — path to save the figure
    """
    valid = [(co2, r) for co2, r in zip(co2_values, results) if r is not None]
    if not valid:
        print("No valid results — skipping phase-evolution plot.")
        return None

    co2_x = np.array([c for c, _ in valid])
    ph_y  = np.array([r['pH'] for _, r in valid])

    # Only show phases with vfrac > 0.5 % somewhere along the CO₂ scan
    all_phases = sorted({
        phase
        for _, r in valid
        for phase, vf in r['phase_vfrac'].items()
        if vf > 0.005
    })

    n_phases = max(len(all_phases), 1)
    colors   = plt.get_cmap('tab20')(np.linspace(0, 1, n_phases))
    hatches  = ['/', '-', '|', '+', 'x', '\\', 'O', '.', '//', '--', '|-', 'XX']

    fig, ax1 = plt.subplots(figsize=(11, 6))
    prev = np.zeros(len(co2_x))
    handles = []

    for i, phase in enumerate(all_phases):
        vals = np.array([r['phase_vfrac'].get(phase, 0.0) for _, r in valid])
        ax1.fill_between(
            co2_x, prev, prev + vals,
            color=colors[i],
            hatch=hatches[i % len(hatches)],
            edgecolor='black', linewidth=0.5,
            label=phase, alpha=0.85,
        )
        prev = prev + vals
        handles.append(mpatches.Patch(
            facecolor=colors[i],
            hatch=hatches[i % len(hatches)],
            edgecolor='black', label=phase,
        ))

    ax1.set_xlabel('Amount of CO$_2$ [g / 100 g cement blend]', fontsize=12)
    ax1.set_ylabel('Volume fraction [m$^3$ / m$^3$ of initial volume]', fontsize=12)
    ax1.set_xlim(CO2_MIN, CO2_MAX)
    ax1.set_ylim(0, 1.15)

    # pH on secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(co2_x, ph_y, 'k-o', linewidth=2, markersize=5, label='pH', zorder=10)
    ax2.set_ylabel('pH', fontsize=12)
    ax2.set_ylim(5, 15)

    ax1.legend(
        handles=handles, ncol=4,
        bbox_to_anchor=(0, 1.02), loc='lower left',
        borderaxespad=0., fontsize='small',
        title='Phases', handlelength=3,
    )
    ax2.legend(loc='upper right', fontsize=10)

    fig.suptitle('Phase Evolution vs CO$_2$ Amount  (Figure 14a)', fontsize=13, y=1.01)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved  →  {save_path}")

    plt.close(fig)
    return fig


def plot_ionic_concentration(co2_values, results, save_path=None):
    """
    Figure 14(b) — semi-log plot of ionic concentrations vs CO₂ amount.

    Parameters
    ----------
    co2_values : array-like
    results    : list of dict (None entries are skipped)
    save_path  : str or None
    """
    valid = [(co2, r) for co2, r in zip(co2_values, results) if r is not None]
    if not valid:
        print("No valid results — skipping ionic-concentration plot.")
        return None

    co2_x      = np.array([c for c, _ in valid])
    line_styles = ['-', '--', '-.', ':', '-', '--', '-.', ':']
    markers     = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
    colors_ion  = plt.get_cmap('tab10')(np.linspace(0, 1, len(IONIC_SPECIES)))

    fig, ax = plt.subplots(figsize=(9, 6))
    plotted = False

    for i, (species, label) in enumerate(IONIC_SPECIES.items()):
        concs = np.array([r['aq_composition'].get(species, 1e-20) for _, r in valid])
        if np.max(concs) < 1e-15:
            continue
        ax.semilogy(
            co2_x, concs,
            linestyle=line_styles[i % len(line_styles)],
            marker=markers[i % len(markers)],
            color=colors_ion[i],
            linewidth=1.8, markersize=5,
            label=label,
        )
        plotted = True

    if not plotted:
        print("No ionic species found with non-negligible concentration.")
        plt.close(fig)
        return None

    ax.set_xlabel('Amount of CO$_2$ [g / 100 g cement blend]', fontsize=12)
    ax.set_ylabel(r'Concentration [mol / kg$_{\mathrm{H_2O}}$]', fontsize=12)
    ax.set_xlim(CO2_MIN, CO2_MAX)
    ax.legend(ncol=2, fontsize=10, loc='best')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_title('Ionic Concentrations vs CO$_2$ Amount  (Figure 14b)', fontsize=13)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved  →  {save_path}")

    plt.close(fig)
    return fig

# ===========================================================================
# Entry point
# ===========================================================================

def main():
    print("=" * 65)
    print("  GEMS CO₂ Carbonation Simulation")
    print("=" * 65)

    total_binder = sum(CLINK_PHASES.values()) + GGBFS_CONTENT
    print("\nFixed cement recipe (g / 100 g cement blend):")
    for phase, amt in CLINK_PHASES.items():
        print(f"  {phase:<8} {amt:6.1f} g")
    print(f"  {'GGBFS':<8} {GGBFS_CONTENT:6.1f} g")
    print(f"  {'Gypsum':<8} {GYPSUM_CONTENT:6.1f} g  (Gp = CaSO₄·2H₂O)")
    print(f"  {'Water':<8} {WATER_CONTENT:6.1f} g  "
          f"(w/b ≈ {WATER_CONTENT / total_binder:.2f})")
    print(f"  Temperature : {TEMPERATURE} °C\n")

    co2_values, results = run_co2_scan()

    n_ok = sum(1 for r in results if r is not None)
    print(f"\n{n_ok} / {len(results)} equilibration steps converged.")

    if n_ok == 0:
        print("No successful equilibrations — aborting.")
        return

    export_csv(co2_values, results)

    print("Generating plots …")
    plot_phase_evolution(co2_values, results,       save_path=FIG_PHASE_OUT)
    plot_ionic_concentration(co2_values, results,   save_path=FIG_IONIC_OUT)

    print("\nAll done.  Output files:")
    print(f"  {CSV_OUTPUT}")
    print(f"  {FIG_PHASE_OUT}")
    print(f"  {FIG_IONIC_OUT}")


if __name__ == "__main__":
    main()
