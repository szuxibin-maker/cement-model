#!/usr/bin/env python3
"""Standalone GEMS CO₂ carbonation simulation script.

Simulates the effect of progressively increasing CO₂ amounts on a fully
hydrated cement paste and generates results similar to Figure 14:

  (a) Two-subplot figure: pH vs CO₂ (top) + stacked solid mass (bottom)
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

# Global font settings (Times New Roman throughout)
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['mathtext.fontset'] = 'stix'

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
    'Al+3':    'Al',
    'CO3-2':   'C',
    'Ca+2':    'Ca',
    'K+':      'K',
    'Mg+2':    'Mg',
    'Na+':     'Na',
    'SO4-2':   'S',
    'H4SiO4@': 'Si',
    'OH-':     'OH-',
}

# Phase merge rules: GEMS phase name → merged display name
MERGE_RULES = {
    # ---- AFt (ettringite family) ------------------------------------------
    "ettringite":       "AFt",
    "C6AsH13":          "AFt",
    "C6AsH9":           "AFt",
    "thaumasite":       "AFt",
    "ettringite-AlFe":  "AFt",
    "ettringite-FeAl":  "AFt",
    "SO4_CO3_AFt":      "AFt",
    "CO3_SO4_AFt":      "AFt",
    # ---- AFm (monosulphate / monocarbonate family) ------------------------
    "C4AH19":           "AFm",
    "C4AH13":           "AFm",
    "C4AH11":           "AFm",
    "CAH10":            "AFm",
    "C4Ac0.5H12":       "AFm",
    "C4Ac0.5H105":      "AFm",
    "C4Ac0.5H9":        "AFm",
    "C4AcH11":          "AFm",
    "C4AcH9":           "AFm",
    "C4AsH16":          "AFm",
    "C4AsH14":          "AFm",
    "C4AsH12":          "AFm",
    "C4AsH105":         "AFm",
    "C4AsH9":           "AFm",
    "Friedels":         "AFm",
    "C2ASH55":          "AFm",
    "C4FH13":           "AFm",
    "C4Fc05H10":        "AFm",
    "C4FcH12":          "AFm",
    "monosulph-FeAl":   "AFm",
    "monosulph-AlFe":   "AFm",
    "SO4_OH_AFm":       "AFm",
    "OH_SO4_AFm":       "AFm",
    "C2AH75":           "AFm",
    "Kuzels":           "AFm",
    # ---- Hydrogarnet -------------------------------------------------------
    "C3AH6":            "Hydrogarnet",
    "C3FH6":            "Hydrogarnet",
    "C3FS0.84H4.32":    "Hydrogarnet",
    "C3(AF)S0.84H":     "Hydrogarnet",
    "C3FS1.34H3.32":    "Hydrogarnet",
    "C3(AF)S0.8":       "Hydrogarnet",
    # ---- Al(OH)3 -----------------------------------------------------------
    "Al(OH)3am":        "Al(OH)3",
    "Al(OH)3mic":       "Al(OH)3",
    "Gibbsite":         "Al(OH)3",
    # ---- Hydrotalcite ------------------------------------------------------
    "hydrotalc-pyro":   "Hydrotalcite",
    "OH-hydrotalcite":  "Hydrotalcite",
    # ---- M-S-H -------------------------------------------------------------
    "MSH":              "M-S-H",
    # ---- Zeolites ----------------------------------------------------------
    "ZeoliteP":         "Zeolites",
    "Natrolite":        "Zeolites",
    "Chabazite":        "Zeolites",
    "ZeoliteX":         "Zeolites",
    "ZeoliteY":         "Zeolites",
    # ---- CSHQ solid solution (species expanded by _solid_masses_g) --------
    "CSHQ-JenD":        "high_Ca_CSH",
    "CSHQ-JenH":        "high_Ca_CSH",
    "CSHQ-TobD":        "high_Ca_CSH",
    "CSHQ-TobH":        "Decalcified_CSH",
    # ---- Straetlingite -----------------------------------------------------
    "straetlingite":    "Straetlingite",
}

# Minimum threshold for displaying a phase (g per 100 g blend)
PHASE_MASS_THRESHOLD_G = 1.0

# ===========================================================================
# Phase suppression  (consistent with util/final_hydration.py)
# ===========================================================================

def _suppress_phases(gemsk):
    """Suppress phases irrelevant for OPC + GGBFS carbonation."""
    _to_suppress = [
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
                'co2':               co2,
                'status':            status,
                'pH':                gemsk.pH,
                'pE':                gemsk.pE,
                'ionic_strength':    gemsk.ionic_strength,
                'phase_volumes':     phase_vols,
                'phase_masses':      gemsk.phase_masses,
                'cshq_species_masses': gemsk.cshq_species_masses,
                'phase_vfrac':       phase_vfrac,
                'aq_composition':    gemsk.aq_composition,
            })
            print(f"OK  (pH = {results[-1]['pH']:.2f})")

        except Exception as exc:
            print(f"ERROR — {exc}")
            results.append(None)

    return co2_values, results

# ===========================================================================
# Helpers
# ===========================================================================

def _solid_masses_g(r):
    """Return a flat dict {phase_or_species_name: mass_g} for all solid phases.

    The aggregate CSHQ solid-solution entry from ``phase_masses`` is replaced
    by the individual CSHQ species masses (``cshq_species_masses``), so that
    ``MERGE_RULES`` can properly split them into *high_Ca_CSH* and
    *Decalcified_CSH*.  All other phases are converted from kg to g.
    """
    masses = {}
    for phase, mass_kg in r['phase_masses'].items():
        if phase in ('aq_gen', 'gas_gen', 'CSHQ'):
            continue
        masses[phase] = mass_kg * 1000.0          # kg → g
    # Expand CSHQ solid solution into individual species (already in grams)
    for species, mass_g in r.get('cshq_species_masses', {}).items():
        masses[species] = mass_g
    return masses


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
        # Phase masses in grams — CSHQ expanded to individual species
        for phase, mass_g in _solid_masses_g(res).items():
            row[f'{phase}_mass'] = mass_g
        # Phase volume fractions
        for phase, vfrac in res['phase_vfrac'].items():
            row[f'vfrac_{phase}'] = vfrac
        # Aqueous concentrations
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
    Figure 14(a) — two-subplot figure:
      top    : pH vs CO₂ amount
      bottom : stacked solid mass (g / 100 g blend) vs CO₂ amount

    Phases are merged via MERGE_RULES, filtered at PHASE_MASS_THRESHOLD_G,
    and drawn with hatches + tab20 colours (Times New Roman font).

    Parameters
    ----------
    co2_values : array-like
    results    : list of dict (None entries are skipped)
    save_path  : str or None
    """
    valid = [(co2, r) for co2, r in zip(co2_values, results) if r is not None]
    if not valid:
        print("No valid results — skipping phase-evolution plot.")
        return None

    co2_x = np.array([c for c, _ in valid])
    ph_y  = np.array([r['pH'] for _, r in valid])

    # ------------------------------------------------------------------
    # Build merged phase mass matrix (grams per 100 g blend)
    # ------------------------------------------------------------------
    phase_series_dict = {}  # {merged_name: np.ndarray of mass values (g)}

    for step_i, (_, r) in enumerate(valid):
        for phase, mass_g in _solid_masses_g(r).items():
            merged = MERGE_RULES.get(phase, phase)
            if merged not in phase_series_dict:
                phase_series_dict[merged] = np.zeros(len(valid))
            phase_series_dict[merged][step_i] += mass_g

    # Filter: skip phases whose maximum mass is below threshold
    phase_series_dict = {
        name: arr
        for name, arr in phase_series_dict.items()
        if arr.max() > PHASE_MASS_THRESHOLD_G
    }

    if not phase_series_dict:
        print("No phases above mass threshold — skipping phase-evolution plot.")
        return None

    # ------------------------------------------------------------------
    # Phase ordering: Decalcified_CSH at index 5, high_Ca_CSH at index 6
    # ------------------------------------------------------------------
    unique_phases = list(phase_series_dict.keys())
    for target in ('Decalcified_CSH', 'high_Ca_CSH'):
        if target in unique_phases:
            unique_phases.remove(target)
    # Re-insert at positions 5 / 6 (clamped to list length)
    for insert_pos, target in ((5, 'Decalcified_CSH'), (6, 'high_Ca_CSH')):
        if target in phase_series_dict:
            unique_phases.insert(min(insert_pos, len(unique_phases)), target)

    # ------------------------------------------------------------------
    # Colours and hatches
    # ------------------------------------------------------------------
    custom_hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*',
                      '//', '--', '||', '++']
    cmap   = plt.get_cmap('tab20')
    colors = cmap(np.linspace(0, 1, max(len(unique_phases), 1)))

    phase_colors  = {p: colors[i]                              for i, p in enumerate(unique_phases)}
    phase_hatches = {p: custom_hatches[i % len(custom_hatches)] for i, p in enumerate(unique_phases)}

    # ------------------------------------------------------------------
    # Two-subplot figure
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        gridspec_kw={'height_ratios': [0.8, 3]},
        figsize=(9, 7),
        sharex=True,
    )

    # ---- Subplot 1: pH vs CO₂ ----------------------------------------
    ax1.plot(co2_x, ph_y, marker='o', linestyle='-', label='pH')
    ax1.set_ylabel('pH', fontsize=20)
    ax1.tick_params(axis='both', labelsize=18)
    ax1.legend(loc='upper right', fontsize=15)
    ax1.grid(True)

    # ---- Subplot 2: stacked solid mass --------------------------------
    prev = np.zeros(len(co2_x))
    for phase in unique_phases:
        vals = phase_series_dict[phase]
        ax2.fill_between(
            co2_x, prev, prev + vals,
            color=phase_colors[phase],
            hatch=phase_hatches[phase],
            edgecolor='black',
            label=phase,
        )
        prev = prev + vals

    # Custom legend patches
    handles = [
        mpatches.Patch(
            facecolor=phase_colors[p],
            hatch=phase_hatches[p],
            edgecolor='black',
            label=p,
        )
        for p in unique_phases
    ]
    ax2.legend(handles=handles, ncol=1,
               bbox_to_anchor=(1, 1), loc='upper left', fontsize=12)

    ax2.set_xlabel('Amount of CO\u2082 (g/100 cement blend)', fontsize=18)
    ax2.set_ylabel('Solid Mass (g)', fontsize=18)
    ax2.tick_params(axis='both', labelsize=18)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved  →  {save_path}")

    plt.close(fig)
    return fig


def plot_ionic_concentration(co2_values, results, save_path=None):
    """
    Figure 14(b) — semi-log plot of ionic concentrations vs CO₂ amount.

    Plots species defined in IONIC_SPECIES with large fonts (Times New Roman)
    and the legend placed outside the axes.

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

    co2_x = np.array([c for c, _ in valid])

    fig, ax = plt.subplots(figsize=(10, 7))
    plotted = False

    for species, label in IONIC_SPECIES.items():
        y = np.array([r['aq_composition'].get(species, 0.0) for _, r in valid])
        mask = y > 0
        if mask.sum() == 0:
            continue
        ax.plot(
            co2_x[mask], y[mask],
            marker='o',
            linestyle='-',
            label=label,
            linewidth=2,
        )
        plotted = True

    if not plotted:
        print("No ionic species found with non-negligible concentration.")
        plt.close(fig)
        return None

    ax.set_xlabel('Amount of CO\u2082 (g/100 cement blend)', fontsize=26)
    ax.set_ylabel('Ions Concentration (M)', fontsize=26)
    ax.tick_params(axis='both', labelsize=24)
    ax.set_yscale('log')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20)
    ax.grid(True, which='major', linestyle='--')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
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
