import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from run.GEMSCalc import GEMS
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cmocean
import os
from run import Time_points

matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['mathtext.fontset'] = 'stix'

# ============================================================
# 相合并规则：把 GEMS 内部多个别名统一到一个显示名
# ============================================================
MERGE_RULES = {
    "SO4_CO3_AFt":    "Ettringite",
    "CO3_SO4_AFt":    "Ettringite",
    "ettringite":     "Ettringite",
    "C3(AF)S0.8":     "Hydrogarnet",
    "C3FS0.84H4.32":  "Hydrogarnet",
    "C3(AF)S0.84H":   "Hydrogarnet",
    "CSHQ-JenD":      "high_Ca_CSH",
    "CSHQ-JenH":      "high_Ca_CSH",
    "CSHQ-TobD":      "high_Ca_CSH",
    "CSHQ-TobH":      "Decalcified_CSH",
    "C4AsH14":        "AFm",
    "straetlingite":  "Straetlingite",
    "C4AsH12":        "AFm",
    "OH_SO4_AFm":     "AFm",
    "C4AcH11":        "Monocarboaluminate",
    "C4Ac0.5H12":     "Hemicarbonate",
    "OH-hydrotalcite": "OH-quintinite",
}

# ============================================================
# 显式物理顺序：底部（最稳定/最主要）→ 顶部（碳化新生）
# ============================================================
PHASE_ORDER = [
    # 未水化熟料
    "C3S",
    "C2S",
    "C3A",
    "C4AF",
    # 主要水化产物 C-S-H
    "high_Ca_CSH",
    "Decalcified_CSH",
    # Portlandite
    "Portlandite",
    # AFt
    "Ettringite",
    # AFm 系列
    "AFm",
    "Monocarboaluminate",
    "Hemicarbonate",
    # 其他水化产物
    "Hydrogarnet",
    "Straetlingite",
    "OH-quintinite",
    # 碳化产物（随 CO2 增加新生，放顶部）
    "Calcite",
    "Silica-amorph",
]

# ============================================================
# 辅助函数：按 PHASE_ORDER 对 unique_phases 排序
# ============================================================
def _sort_phases(phases_iterable):
    """
    按 PHASE_ORDER 排序，不在列表中的相追加到末尾（保持稳定顺序）。
    """
    phases = list(phases_iterable)
    ordered = [p for p in PHASE_ORDER if p in phases]
    extra   = [p for p in phases if p not in PHASE_ORDER]
    return ordered + extra


# ============================================================
# Parrot & Killoh 模型参数
# ============================================================
K1 = {"C3S":1.5,"C2S":0.5,"C3A":1.,"C4AF":0.37}
N1 = {"C3S":0.7,"C2S":1.0,"C3A":0.85,"C4AF":0.7}
K2 = {"C3S":0.05,"C2S":0.02,"C3A":0.04,"C4AF":0.015}
K3 = {"C3S":1.1,"C2S":0.7,"C3A":1.0,"C4AF":0.4}
N3 = {"C3S":3.3,"C2S":5.0,"C3A":3.2,"C4AF":3.7}
H  = {"C3S":1.8,"C2S":1.35,"C3A":1.6,"C4AF":1.45}
Ea = {"C3S":41570,"C2S":20785,"C3A":54040,"C4AF":34087}
T0 = 25  # °C
ref_fineness = 385

output_times = Time_points.output_times


def bouges_composition(CaO, SiO2, Al2O3, Fe2O3, SO3):
    """
    Computation of clinker phases from oxides using Bogue's conversion.
    """
    clink_phases = {}
    clink_phases["C3S"]  = C3S  = 4.07*CaO - 7.6*SiO2 - 6.72*Al2O3 - 1.43*Fe2O3 - 2.85*SO3
    clink_phases["C2S"]  = C2S  = 2.87*SiO2 - 0.75*C3S
    clink_phases["C3A"]  = C3A  = 2.65*Al2O3 - 1.69*Fe2O3
    clink_phases["C4AF"] = C4AF = 3.04*Fe2O3
    return clink_phases


"""
Parrot & Killoh rate model
"""
def nuclandgrowth(alpha_t, K1, N1, fineness, ref_fineness):
    ln = np.log
    rt = (K1/N1) * (1-alpha_t) * (-ln(1-alpha_t))**(1-N1) * (fineness/ref_fineness)
    return rt

def diffusion(alpha_t, K2):
    eps = 1e-12
    alpha_safe = np.clip(alpha_t, eps, 1 - eps)
    rt = K2 * (1 - alpha_safe)**(2/3) / (alpha_safe)**(1/3)
    return rt

def hydrationshell(alpha_t, K3, N3):
    rt = K3 * (1-alpha_t)**N3
    return rt

def f_wc_lothenbach(wc, H, alpha_t):
    alpha_cr = H * wc
    if alpha_t > alpha_cr:
        return (1 + 3.333*(H*wc - alpha_t))**4
    else:
        return 1

def f_rh(RH):
    return ((RH - 0.55)/0.45)**4

def f_temp(T, T0, Ea):
    exp = np.exp
    R   = 8.314
    T0  = T + 273.15
    T   = T + 273.15
    return exp((-Ea/R) * ((1/T) - (1/T0)))

def overall_rate(t, alpha_t, K1, N1, K2, K3, N3, H, wc, RH, fineness, ref_fineness, T, T0, Ea):
    if alpha_t >= 1:
        alpha_t = 0.9999
    r1 = nuclandgrowth(alpha_t, K1, N1, fineness, ref_fineness)
    r2 = diffusion(alpha_t, K2)
    r3 = hydrationshell(alpha_t, K3, N3)
    r  = min(r1, r2, r3)
    r  = r * f_wc_lothenbach(wc, H, alpha_t) * f_rh(RH) * f_temp(T, T0, Ea)
    return r

def parrot_killoh(wc, RH, T, fineness):
    DoH = {}
    for phase in ["C3S","C2S","C3A","C4AF"]:
        ivp_out = solve_ivp(
            overall_rate, [0, np.max(output_times)], [1e-15],
            t_eval=output_times,
            args=(K1[phase], N1[phase], K2[phase], K3[phase], N3[phase],
                  H[phase], wc, RH, fineness, ref_fineness, T, T0, Ea[phase])
        )
        DoH[phase] = np.squeeze(ivp_out.y[0])
    return DoH


def plot_bars(results):
    """Bar chart of phase volumes/masses at each time step."""
    unique_phases = []
    for t in results.keys():
        for phase in results[t].keys():
            if results[t][phase] > 0:
                if phase not in unique_phases:
                    unique_phases.append(phase)
    colors = plt.get_cmap('rainbow')(np.linspace(0., 1, len(unique_phases)+1))
    phase_colors = {}
    for phase, color in zip(unique_phases, colors):
        phase_colors[phase] = color
    i = 0
    ytks = []
    ytksval = []
    for t in results.keys():
        prev = 0
        for key, val in results[t].items():
            if key in unique_phases:
                plt.barh(i, val, left=prev, height=0.5,
                         label=key, color=phase_colors[key])
                prev += val
        ytks.append(i)
        ytksval.append(str(t))
        i += 1
        if t == 0:
            plt.legend(ncol=3, bbox_to_anchor=(0, 1.5), loc='upper left')
    plt.yticks(ytks, ytksval)
    plt.xlabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    plt.ylabel("Time [days]")
    return plt


def to_phase_first_dict(results):
    """Convert time-first dict to phase-first dict for phase diagrams."""
    out = {}
    for t in results.keys():
        for phase in results[t].keys():
            if phase not in out:
                out[phase] = []
            out[phase].append(results[t][phase])
    for phase in out.keys():
        out[phase] = np.array(out[phase])
    return out


def phase_plot(results, datatype="vfrac"):
    """
    Plot phase diagrams vs. time with controlled phase order and improved aesthetics.

    Parameters
    ----------
    results  : dict  {phase: np.array(values_over_time)}
    datatype : str   "vfrac" or "masses"
    """
    # --- 筛选有效相 ---
    raw_phases = []
    for phase in results.keys():
        if np.max(results[phase]) > 0.01:
            raw_phases.append(phase)

    # --- 按物理顺序排列 ---
    unique_phases = _sort_phases(raw_phases)

    # --- 颜色与填充图案 ---
    custom_hatches = ['/', '\', '|', '-', '+', 'x', 'o', 'O', '.', '*',
                      '//', '\\', '||', '--', '++']
    cmap   = plt.get_cmap('tab20')
    colors = cmap(np.linspace(0, 1, max(len(unique_phases), 1)))
    phase_hatches = {p: custom_hatches[i % len(custom_hatches)]
                     for i, p in enumerate(unique_phases)}
    phase_colors  = {p: colors[i] for i, p in enumerate(unique_phases)}

    # --- 绘图 ---
    fig, ax = plt.subplots(figsize=(9, 6))

    prev = np.zeros(len(output_times))
    for phase in unique_phases:
        vals = results[phase]
        ax.fill_between(output_times, prev, prev + vals,
                        color=phase_colors[phase],
                        hatch=phase_hatches[phase],
                        edgecolor="black", linewidth=0.5,
                        label=phase)
        prev = prev + vals

    handles = [
        mpatches.Patch(facecolor=phase_colors[p], hatch=phase_hatches[p],
                       edgecolor='black', label=p)
        for p in unique_phases
    ]
    ax.legend(handles=handles, ncol=1,
              bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=11, title='Phases', title_fontsize=12)

    if datatype == "vfrac":
        ax.set_ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]', fontsize=16)
        ax.set_ylim(0, 1)
    elif datatype == "masses":
        ax.set_ylabel('Solid Mass (g)', fontsize=16)

    ax.set_xlabel("Time [days]", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    return fig, ax


def co2_phase_plot(results, datatype="masses",
                   save_path=None, target_ph=None, co2_at_target=None):
    """
    Plot stacked phase diagram and optional pH curve vs. CO2 amount.

    Parameters
    ----------
    results        : dict  {co2_val: {phase: scalar}}
                     OR a pandas DataFrame with columns [CO2, pH, *_mass/*_vfrac]
    datatype       : str   "masses" or "vfrac"
    save_path      : str   optional output file path (e.g. './picture results/out.png')
    target_ph      : float optional target pH line drawn on pH subplot
    co2_at_target  : float optional CO2 value annotation at target pH

    Notes
    -----
    When `results` is a DataFrame it must contain:
      - 'CO2' column
      - 'pH'  column  (for the top sub-plot)
      - columns ending with '_mass' (datatype="masses") or '_vfrac' (datatype="vfrac")
    The MERGE_RULES dict is applied to consolidate phase aliases, and phases are
    rendered in the canonical PHASE_ORDER (bottom → top).
    """

    # ------------------------------------------------------------------ #
    # 1.  统一数据格式：把 dict 和 DataFrame 都转成
    #     phase_series_dict  {merged_phase: pd.Series(index=row_index)}
    #     co2_values         sorted array
    #     ph_series          pd.Series or None
    # ------------------------------------------------------------------ #
    if isinstance(results, pd.DataFrame):
        df = results.copy()
        co2_values = np.sort(df['CO2'].values)
        df = df.set_index('CO2').sort_index()
        ph_series = df['pH'] if 'pH' in df.columns else None

        suffix = '_mass' if datatype == "masses" else '_vfrac'
        threshold = 1 if datatype == "masses" else 0.01

        phase_series_dict = {}
        for col in df.columns:
            if not col.endswith(suffix):
                continue
            raw = col[:-len(suffix)]
            if raw in ("aq_gen", "gas_gen"):
                continue
            merged = MERGE_RULES.get(raw, raw)
            if df[col].max() <= threshold:
                continue
            if merged not in phase_series_dict:
                phase_series_dict[merged] = df[col].copy()
            else:
                phase_series_dict[merged] = phase_series_dict[merged] + df[col]

        # build value array per phase aligned to co2_values
        def _get_vals(phase):
            s = phase_series_dict[phase]
            return np.array([s.get(c, 0.0) for c in co2_values])

    else:
        # legacy dict path: {co2_val: {phase: scalar}}
        co2_values = np.array(sorted(results.keys()))
        ph_series  = None
        threshold  = 1 if datatype == "masses" else 0.01

        raw_phase_set = {}
        for co2, data in results.items():
            for phase, val in data.items():
                if val > threshold:
                    merged = MERGE_RULES.get(phase, phase)
                    if merged not in raw_phase_set:
                        raw_phase_set[merged] = {}
                    raw_phase_set[merged][co2] = raw_phase_set[merged].get(co2, 0) + val

        phase_series_dict = raw_phase_set  # {merged: {co2: val}}

        def _get_vals(phase):
            d = phase_series_dict[phase]
            return np.array([d.get(c, 0.0) for c in co2_values])

    if not phase_series_dict:
        print("没有满足条件的相可以绘图。")
        return None, None

    # ------------------------------------------------------------------ #
    # 2.  按 PHASE_ORDER 排序
    # ------------------------------------------------------------------ #
    unique_phases = _sort_phases(phase_series_dict.keys())

    # ------------------------------------------------------------------ #
    # 3.  颜色与填充图案
    # ------------------------------------------------------------------ #
    custom_hatches = ['/', '\', '|', '-', '+', 'x', 'o', 'O', '.', '*',
                      '//', '\\', '||', '--', '++']
    cmap   = plt.get_cmap('tab20')
    colors = cmap(np.linspace(0, 1, max(len(unique_phases), 1)))
    phase_hatches = {p: custom_hatches[i % len(custom_hatches)]
                     for i, p in enumerate(unique_phases)}
    phase_colors  = {p: colors[i] for i, p in enumerate(unique_phases)}

    # ------------------------------------------------------------------ #
    # 4.  建图：有 pH 时用双子图，否则单图
    # ------------------------------------------------------------------ #
    if ph_series is not None:
        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            gridspec_kw={'height_ratios': [0.8, 3]},
            figsize=(9, 7),
            sharex=True
        )
        # --- 子图 1: pH ---
        ph_co2   = ph_series.index.values
        ph_vals  = ph_series.values
        ax1.plot(ph_co2, ph_vals, marker='o', linestyle='-',
                 linewidth=2, label='pH vs CO₂')
        ax1.set_ylabel('pH', fontsize=20)
        ax1.tick_params(axis='both', labelsize=18)
        ax1.grid(True, linestyle='--', alpha=0.4)
        if target_ph is not None:
            ax1.axhline(y=target_ph, linestyle='--', color='green',
                        linewidth=1.5, label=f'pH = {target_ph}')
        ax1.legend(loc='upper right', fontsize=14)
        ax_stack = ax2
    else:
        fig, ax_stack = plt.subplots(figsize=(9, 6))

    # --- 子图（或主图）: 堆叠面积 ---
    prev = np.zeros(len(co2_values))
    for phase in unique_phases:
        vals = _get_vals(phase)
        ax_stack.fill_between(
            co2_values, prev, prev + vals,
            color=phase_colors[phase],
            hatch=phase_hatches[phase],
            edgecolor="black",
            label=phase
        )
        prev = prev + vals

    handles = [
        mpatches.Patch(facecolor=phase_colors[p], hatch=phase_hatches[p],
                       edgecolor='black', label=p)
        for p in unique_phases
    ]
    ax_stack.legend(handles=handles, ncol=1,
                    bbox_to_anchor=(1.01, 1), loc='upper left',
                    fontsize=12, title='Phases', title_fontsize=13)

    ax_stack.set_xlabel("Amount of CO₂ (g/100 g cement blend)", fontsize=18)
    if datatype == "masses":
        ax_stack.set_ylabel('Solid Mass (g)', fontsize=18)
    else:
        ax_stack.set_ylabel('Volume fraction [m³/m³]', fontsize=18)
    ax_stack.tick_params(axis='both', labelsize=16)
    ax_stack.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图已保存到：{save_path}")

    plt.show()
    return fig, ax_stack


def run_hydration(clink_phases, wc, CSH2, T, DoH):
   input_file = 'gems_files/CemHyds-dat.lst'
   gemsk = GEMS(input_file)
   all_species = clink_phases.copy()
   all_species["H2O@"] = wc * 100
   all_species["Gp"] = CSH2
   for name in all_species.keys():
       all_species[name] *= 1e-3
   gemsk.T = T + 273.15
   gemsk.add_multiple_species_amt(all_species, units="kg")
   gemsk.add_species_amt("O2", 1e-6)
   gems_vol_frac = {}
   gems_phase_masses = {}
   density = []
   for i in range(len(output_times)):
       if i > 0: gemsk.warm_start()
       for phase in clink_phases:
           gemsk.species_lower_bound(phase, clink_phases[phase]*(1-DoH[phase][i])*1e-3, units="kg")
       print("Time-->", str(output_times[i]), ":", gemsk.equilibrate())
       gems_phase_masses[output_times[i]] = gemsk.phase_masses
       if output_times[i] == 0: init_vol = gemsk.system_volume
       print(init_vol)
       gems_vol_frac[output_times[i]] = gemsk.phase_volumes
       for key, val in gems_vol_frac[output_times[i]].items():
           gems_vol_frac[output_times[i]][key] /= init_vol
       density.append(gemsk.system_mass / gemsk.system_volume)
   return gems_vol_frac, gems_phase_masses, density