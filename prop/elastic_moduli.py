"""
两阶段 Mori–Tanaka（DS）均化

阶段1：CSH + 水化产物 + 水(aq_gen) [+ 可选孔隙Pore] → 起点基体（可选软化：'voigt'/'reuss'/'hill'/'mix'，mix用soft_bias调节软化程度）
阶段2：未水化相 → 按等效杨氏模量从低到高依次加入（微分体积分）

用法示例：
df_out = add_moduli_two_stage(
    df,
    start_mode='mix',   # 'voigt' | 'reuss' | 'hill' | 'mix'
    soft_bias=0.8,      # 仅在 mode='mix' 时生效，越大越偏软
    include_pore_in_stage1=True,  # 阶段1是否把Pore纳入基体
    n_steps_per_inclusion=50
)
"""
import warnings
warnings.filterwarnings("ignore")


import pandas as pd
import numpy as np

# ===== 模量表 =====
def build_moduli_tables_no_time():
    bulk = {
        # —— 水化产物（不含CSH与水）——
        'SO4_CO3_AFt': 26, 'CO3_SO4_AFt': 26, 'ettringite-AlFe': 26,
        'ettringite-FeAl': 26, 'ettringite': 26, 'straetlingite': 23,
        'SO4_OH_Afm': 25.2, 'OH_SO4_AFm': 25.2, 'C4AsH14': 25.2,
        'C4Ac0.5H12': 25.2, 'C4AcH11': 40, 'Portlandite': 40,
        'C3AH6': 99.8, 'C3(AF)S0.84H': 99.8, 'OH-hydrotalcite': 25.2,
        'C3FS1.34H3.32': 99.8, 'MSH': 14.9,
        # 一些可能归为“水化/碳化产物或填料”的相
        'limestone': 69.8, 'Calcite': 69.8, 'Silica-amorph': 31.2,
        'calcined_clay': 24.39,  # 若你认为它未水化，可在下面集合调整
        'Gypsum': 44.8, 'Gp': 44.8,
        # —— 水、CSH、孔隙 —— 
        'aq_gen': 2,
        'LD_CSH': 13.9, 'HD_CSH': 18.8,
        'Pore': 0,
        # —— 未水化相（熟料/掺合料等）——
        'GGBFS': 38.0, 'fly_ash': 58.8, 'silica_fume': 36.5,
        'Alite': 105.2, 'Belite': 105.2, 'Aluminate': 105.2, 'Ferrite': 105.2,
        # —— 兜底 —— 
        'other_material': 13.9,
    }
    shear = {
        'SO4_CO3_AFt': 9.3, 'CO3_SO4_AFt': 9.3, 'ettringite-AlFe': 9.3,
        'ettringite-FeAl': 9.3, 'ettringite': 9.3, 'straetlingite': 10.6,
        'SO4_OH_Afm': 11.6, 'OH_SO4_AFm': 11.6, 'C4AsH14': 11.6,
        'C4Ac0.5H12': 11.6, 'C4AcH11': 16.0, 'Portlandite': 16.0,
        'C3AH6': 64.3, 'C3(AF)S0.84H': 64.3, 'OH-hydrotalcite': 11.6,
        'C3FS1.34H3.32': 64.3, 'MSH': 9.0,
        'limestone': 30.4, 'Calcite': 30.4, 'Silica-amorph': 36.5,
        'calcined_clay': 18.29,
        'Gypsum': 17.2, 'Gp': 17.2,
        'LD_CSH': 8.8, 'HD_CSH': 11.8, 'aq_gen': 0, 'Pore': 0.0,
        'GGBFS': 28.5, 'fly_ash': 27.2, 'silica_fume': 31.2,
        'Alite': 44.8, 'Belite': 44.8, 'Aluminate': 44.8, 'Ferrite': 44.8,
        'other_material': 8.8,
    }
    return bulk, shear

# ===== 单步 Mori–Tanaka =====
def _mt_update(m_k, m_g, i_k, i_g, i_v):
    if i_v <= 0:
        return m_k, m_g
    denom_m = 3*m_k + m_g
    if denom_m <= 0 or m_k <= 0 or m_g <= 0:
        return m_k, m_g
    nu_m = (3*m_k - 2*m_g) / (2*denom_m)
    alpha = (1 + nu_m) / (3 * (1 - nu_m))
    beta  = 2 * (4 - 5*nu_m) / (15 * (1 - nu_m))
    k_eff = m_k + i_v * (i_k - m_k) / (1 + alpha * (i_k/m_k - 1))
    g_eff = m_g + i_v * (i_g - m_g) / (1 + beta  * (i_g/m_g - 1))
    return k_eff, g_eff

# ===== 组平均（支持软化模式） =====
def _group_avg_KG_vf(vol_phase, group_set, bulk, shear,
                     fallback='other_material', mode='mix', soft_bias=0.7):
    """
    mode: 'voigt'（偏硬）| 'reuss'（偏软）| 'hill'（折中）| 'mix'（Voigt-Reuss 插值）
    soft_bias: 仅在 'mix' 时生效，0→Voigt，1→Reuss
    """
    v = 0.0
    k_sum = 0.0; g_sum = 0.0
    k_inv_sum = 0.0; g_inv_sum = 0.0
    eps = 1e-12
    for ph, vf in vol_phase.items():
        if vf <= 0 or ph not in group_set:
            continue
        k = max(bulk.get(ph, bulk[fallback]), eps)
        g = max(shear.get(ph, shear[fallback]), eps)
        v += vf
        k_sum += k * vf; g_sum += g * vf
        k_inv_sum += vf / k; g_inv_sum += vf / g
    if v <= 0:
        return 0.0, 0.0, 0.0

    k_voigt = k_sum / v
    g_voigt = g_sum / v
    k_reuss = v / k_inv_sum
    g_reuss = v / g_inv_sum

    if mode == 'voigt':
        return k_voigt, g_voigt, v
    elif mode == 'reuss':
        return k_reuss, g_reuss, v
    elif mode == 'hill':
        return 0.5*(k_voigt+k_reuss), 0.5*(g_voigt+g_reuss), v
    elif mode == 'mix':
        return ((1-soft_bias)*k_voigt + soft_bias*k_reuss,
                (1-soft_bias)*g_voigt + soft_bias*g_reuss, v)
    else:
        raise ValueError("Unknown mode")

def _young_from_KG(K, G):
    denom = 3*K + G
    return 0.0 if denom <= 0 else 9*K*G/denom

# ===== 两阶段 MT =====
def compute_two_stage_MT(
    vol_phase: dict, bulk: dict, shear: dict,
    include_pore_in_stage1: bool = True,
    start_mode: str = 'mix', soft_bias: float = 0.7,
    n_steps_per_inclusion: int = 50
):
    """
    阶段1：CSH + 水化产物 + 水（可选含Pore） → 起点基体（可软化）
    阶段2：未水化相 → 按 E 从低到高依次加入
    """
    # ---- 定义两阶段相集合 ----
    STAGE1_SET = {
        # CSH & 水
        'LD_CSH','HD_CSH','aq_gen','Pore'
        # 典型水化/碳铝酸盐/层状产物 & 一些填料类你可按需调整
        # 'Portlandite','ettringite','ettringite-AlFe','ettringite-FeAl',
        # 'SO4_CO3_AFt','CO3_SO4_AFt','straetlingite',
        # 'SO4_OH_Afm','OH_SO4_AFm','C4AsH14','C4Ac0.5H12','C4AcH11',
        # 'C3AH6','C3(AF)S0.84H','OH-hydrotalcite','C3FS1.34H3.32','MSH',
       
        # 视体系而定：可按“是否参与早期水化/碳化沉积/填充”决定是否纳入阶段1
     # 'Silica-amorph',
    }
    if include_pore_in_stage1:
        STAGE1_SET.add('Pore')

    STAGE2_SET = {
        # 未水化熟料 + 常见未反应掺合料
        'Alite','Belite','Aluminate','Ferrite',  'limestone','Calcite','Silica-amorph','calcined_clay',
        'fly_ash','GGBFS','silica_fume','Gypsum',
        'Portlandite','ettringite','ettringite-AlFe','ettringite-FeAl',
        'SO4_CO3_AFt','CO3_SO4_AFt','straetlingite',
        'SO4_OH_Afm','OH_SO4_AFm','C4AsH14','C4Ac0.5H12','C4AcH11',
        'C3AH6','C3(AF)S0.84H','OH-hydrotalcite','C3FS1.34H3.32','MSH',
        # 如果你将 'calcined_clay' 视为未水化，也可以把它移到这里
    }

    # ---- 阶段1：起点（可软化）----
    k1, g1, v1 = _group_avg_KG_vf(vol_phase, STAGE1_SET, bulk, shear,
                                  mode=start_mode, soft_bias=soft_bias)
    if v1 <= 0:
        # 若阶段1集合体积分为0，则无法建立起点
        return 0.0, 0.0
    m_k, m_g = k1, g1
    # print(f"[阶段1起点] K={m_k:.4f}, G={m_g:.4f}, E={_young_from_KG(m_k,m_g):.4f} GPa "
    #       f"(mode={start_mode}, soft_bias={soft_bias})")

    # ---- 阶段2：未水化相按 E 从低到高加入 ----
    inclusions = []
    for ph in STAGE2_SET:
        vf = vol_phase.get(ph, 0.0)
        if vf <= 0:
            continue
        k_i = bulk.get(ph, bulk['other_material'])
        g_i = shear.get(ph, shear['other_material'])
        denom = 3*k_i + g_i
        E_i = 0.0 if denom <= 0 else 9*k_i*g_i/denom
        inclusions.append((E_i, ph, vf, k_i, g_i))
    inclusions.sort(key=lambda x: x[0], reverse=True) # 软 -> 硬

    for E_i, ph, vf, k_i, g_i in inclusions:
        steps = max(int(n_steps_per_inclusion), 1)
        dv = vf / steps
        for _ in range(steps):
            m_k, m_g = _mt_update(m_k, m_g, k_i, g_i, dv)
        # print(f"[阶段2加入] {ph:>15s} vf={vf:.4f} → E={_young_from_KG(m_k,m_g):.4f} GPa")

    return m_k, m_g

# ===== 外层：按 DataFrame 逐行计算 =====
def data_elastic_moduli(
    df_full: pd.DataFrame,
    start_mode: str = 'mix', soft_bias: float = 0.7,
    include_pore_in_stage1: bool = True,
    n_steps_per_inclusion: int = 50
) -> pd.DataFrame:
    """
    df_full：包含若干 `*_volume_frac` 列；可包含 'total_volume_frac_ref0' 用于补干孔
    返回：在副本上新增 K_bulk, G_shear, E_Young, nu_Poisson
    """
    df = df_full.copy()
    bulk, shear = build_moduli_tables_no_time()

    # 仅接受 *_volume_frac，排除总和/气体/已有孔
    def is_frac_candidate(col: str) -> bool:
        if not col.endswith('_volume_frac'):
            return False
        bad = {'total_volume_frac_ref0', 'gas_gen_volume_frac',
               'Pore_volume_frac', 'Vgp_volume_frac', 'Vcp_volume_frac'}
        return col not in bad

    K_list, G_list, E_list, nu_list = [], [], [], []

    for _, r in df.iterrows():
        # 1) 收集各相体积分数
        vol_phase = {}
        for col in df.columns:
            if is_frac_candidate(col):
                f = float(r.get(col, 0.0))
                if f > 0:
                    ph = col[:-12]  # 去掉 "_volume_frac"
                    vol_phase[ph] = vol_phase.get(ph, 0.0) + f

        # 2) 自动补“干孔”：Pore = 1 - total_volume_frac_ref0
        total_ref0 = float(r.get('total_volume_frac_ref0', sum(vol_phase.values())))
        pore = max(0.0, min(1.0, 1.0 - total_ref0))
        if pore > 0:
            vol_phase['Pore'] = vol_phase.get('Pore', 0.0) + pore

        # 3) 两阶段均化
        K, G = compute_two_stage_MT(
            vol_phase, bulk, shear,
            include_pore_in_stage1=include_pore_in_stage1,
            start_mode=start_mode, soft_bias=soft_bias,
            n_steps_per_inclusion=n_steps_per_inclusion
        )

        denom = 3*K + G
        if denom <= 0:
            E = 0.0; nu = 0.0
        else:
            E = 9*K*G/denom
            nu = (3*K - 2*G)/(2*denom)

        K_list.append(K); G_list.append(G); E_list.append(E); nu_list.append(nu)

    df['K_bulk'] = K_list
    df['G_shear'] = G_list
    df['E_Young'] = E_list
    df['nu_Poisson'] = nu_list
    return df
