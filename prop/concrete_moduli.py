import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# # 将弹性模量和泊松比转换为体积模量和剪切模量
# def E_v_to_K_G(E, v):
#     K = E / (3 * (1 - 2 * v))  # 体积模量 K
#     G = E / (2 * (1 + v))      # 剪切模量 G
#     return K, G

# # 将体积模量和剪切模量转换为弹性模量和泊松比
# def K_G_to_E_v(K, G):
#     E = 9 * K * G / (3 * K + G)  # 弹性模量
#     v = (3 * K - 2 * G) / (2 * (3 * K + G))  # 泊松比
#     return E, v

# # 自洽方法计算函数
# def self_consistent_composite(mortar_E, mortar_v, elastic_properties, vol_fractions):
#     """
#     计算基于自洽方法的复合材料有效模量。
    
#     参数:
#     - mortar_E: 基体 (mortar) 的弹性模量 (GPa)
#     - mortar_v: 基体 (mortar) 的泊松比
#     - elastic_properties: 各相的弹性模量和泊松比的字典
#     - vol_fractions: 各相的体积分数字典
    
#     返回:
#     - K_eff: 整体体积模量 (GPa)
#     - G_eff: 整体剪切模量 (GPa)
#     - E_eff: 整体弹性模量 (GPa)
#     - v_eff: 整体泊松比
#     """
    
#     # 初始化基体为 mortar，计算基体的体积模量和剪切模量
#     K_m, G_m = E_v_to_K_G(mortar_E, mortar_v)
    
#     # 自洽方程
#     def equations(x):
#         K_eff, G_eff = x
#         nominator_bulk = 0
#         denominator_bulk = 0
#         nominator_shear = 0
#         denominator_shear = 0
        
#         # 遍历每个相，计算自洽条件下的体积模量和剪切模量
#         for phase, vol_frac in vol_fractions.items():
#             K_i, G_i = E_v_to_K_G(elastic_properties[phase]['E'], elastic_properties[phase]['v'])
            
#             # 自洽方程中的体积模量部分
#             alpha_bulk = (3 * K_eff) / (3 * K_eff + 4 * G_eff)
#             nominator_bulk += vol_frac * K_i * (1 + alpha_bulk * (K_i / K_eff - 1))
#             denominator_bulk += vol_frac * (1 + alpha_bulk * (K_i / K_eff - 1))
            
#             # 自洽方程中的剪切模量部分
#             alpha_shear = (6 * (K_eff + 2 * G_eff)) / (5 * (3 * K_eff + 4 * G_eff))
#             nominator_shear += vol_frac * G_i * (1 + alpha_shear * (G_i / G_eff - 1))
#             denominator_shear += vol_frac * (1 + alpha_shear * (G_i / G_eff - 1))

#         return [K_eff - nominator_bulk / denominator_bulk, G_eff - nominator_shear / denominator_shear]
    
#     # 使用初始猜测的体积模量和剪切模量来求解自洽方程
#     initial_guess = [K_m, G_m]
#     K_eff, G_eff = fsolve(equations, initial_guess)
    
#     # 通过体积模量和剪切模量计算弹性模量和泊松比
#     E_eff, v_eff = K_G_to_E_v(K_eff, G_eff)

#     return K_eff, G_eff, E_eff, v_eff

# # 计算复合材料整体的体积模量、剪切模量、弹性模量和泊松比
# def calculate_composite_moduli(mortar_E, mortar_v):
#     """
#     计算复合材料整体的体积模量、剪切模量、弹性模量和泊松比。
#     """
    
#     # 定义不同相的弹性参数 (E: 弹性模量, v: 泊松比)
#     elastic_properties = {
#         'Stone': {'E': 51.31, 'v': 0.15},   # 骨料
#        # 'ITZ-1': {'E': 10.17, 'v': 0.30},   # ITZ-1 (骨料和砂子之间)
#         'ITZ-2': {'E': 12.7, 'v': 0.30},    # ITZ-2 (砂浆和砂子之间)
#         'Mortar': {'E': mortar_E, 'v': mortar_v}  # mortar
#     }

#     # 定义不同相的体积分数
#     vol_fractions = {
#         'Stone': 0.403,
#         'ITZ-2': 0.151,
#         'Mortar': 0.423
#     }

#     # 使用自洽方法计算复合材料的有效模量
#     K_hom, G_hom, E_hom, v_hom = self_consistent_composite(mortar_E, mortar_v, elastic_properties, vol_fractions)
    
#     return K_hom, G_hom, E_hom, v_hom

# # 处理 pandas 列表中的数据
# def calculate_concrete_moduli(df):
#     """
#     计算 DataFrame 中每一行的体积模量、剪切模量、弹性模量和泊松比，并返回包含计算结果的 DataFrame。
    
#     参数：
#     - df: 包含 mortar_E 和 mortar_v 的 pandas DataFrame。
    
#     返回：
#     - 包含 K_bulk_modulus, G_shear_modulus, E_Young_modulus, v_Poisson_ratio 的 DataFrame。
#     """
#     results = []
#     for index, row in df.iterrows():
#         mortar_E = row['mortar_E']
#         mortar_v = row['mortar_v']
        
#         # 调用计算函数
#         K_hom, G_hom, E_hom, v_hom = calculate_composite_moduli(mortar_E, mortar_v)
        
#         # 保存结果
#         results.append([K_hom, G_hom, E_hom, v_hom])
    
#     # 返回包含结果的 DataFrame
#     return pd.DataFrame(results, columns=['concrete_K', 'concrete_G', 'concrete_E', 'concrete_v'])
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# 将弹性模量和泊松比转换为体积模量和剪切模量
def E_v_to_K_G(E, v):
    K = E / (3 * (1 - 2 * v))  # 体积模量 K
    G = E / (2 * (1 + v))      # 剪切模量 G
    return K, G

# 将体积模量和剪切模量转换为弹性模量和泊松比
def K_G_to_E_v(K, G):
    E = 9 * K * G / (3 * K + G)  # 弹性模量
    v = (3 * K - 2 * G) / (2 * (3 * K + G))  # 泊松比
    return E, v

# Mori-Tanaka 方法计算函数
def mori_tanaka_composite(mortar_E, mortar_v, elastic_properties, vol_fractions):
    """
    计算基于 Mori-Tanaka 方法的复合材料有效模量。
    
    参数:
    - mortar_E: 基体 (mortar) 的弹性模量 (GPa)
    - mortar_v: 基体 (mortar) 的泊松比
    - elastic_properties: 各相的弹性模量和泊松比的字典
    - vol_fractions: 各相的体积分数字典
    
    返回:
    - K_eff: 整体体积模量 (GPa)
    - G_eff: 整体剪切模量 (GPa)
    - E_eff: 整体弹性模量 (GPa)
    - v_eff: 整体泊松比
    """
    
    # 计算基体 (mortar) 的体积模量和剪切模量
    K_m, G_m = E_v_to_K_G(mortar_E, mortar_v)

    # Mori-Tanaka 方法核心计算
    def mori_tanaka(K_m, G_m, K_i, G_i, f_i):
        # 计算 Mori-Tanaka 张量
        S_bulk = (3 * K_m) / (3 * K_m + 4 * G_m)  # 体积模量的 Eshelby 张量
        S_shear = (6 * (K_m + 2 * G_m)) / (5 * (3 * K_m + 4 * G_m))  # 剪切模量的 Eshelby 张量

        # Mori-Tanaka 体积模量
        K_eff = K_m + f_i * (K_i - K_m) * S_bulk / (1 + f_i * S_bulk * (K_i / K_m - 1))
        
        # Mori-Tanaka 剪切模量
        G_eff = G_m + f_i * (G_i - G_m) * S_shear / (1 + f_i * S_shear * (G_i / G_m - 1))
        
        return K_eff, G_eff

    # 初始有效模量为基体模量
    K_eff = K_m
    G_eff = G_m

    # 遍历夹杂相，使用 Mori-Tanaka 方法计算复合材料模量
    for phase, vol_frac in vol_fractions.items():
        if phase == 'Mortar':  # 跳过基体
            continue
        
        # 获取夹杂物的弹性模量和泊松比，并转换为体积模量和剪切模量
        E_i = elastic_properties[phase]['E']
        v_i = elastic_properties[phase]['v']
        K_i, G_i = E_v_to_K_G(E_i, v_i)

        # 通过 Mori-Tanaka 方法更新基体模量
        K_eff, G_eff = mori_tanaka(K_eff, G_eff, K_i, G_i, vol_frac)

    # 通过体积模量和剪切模量计算弹性模量和泊松比
    E_eff, v_eff = K_G_to_E_v(K_eff, G_eff)

    return K_eff, G_eff, E_eff, v_eff

# 计算复合材料整体的体积模量、剪切模量、弹性模量和泊松比
def calculate_composite_moduli(mortar_E, mortar_v):
    """
    计算复合材料整体的体积模量、剪切模量、弹性模量和泊松比。
    """
    
    # 定义不同相的弹性参数 (E: 弹性模量, v: 泊松比)
    elastic_properties = {
        'Stone': {'E': 51.31, 'v': 0.15},   # 骨料
        'ITZ-2': {'E': 12.7, 'v': 0.30},    # ITZ-2 (砂浆和砂子之间)
        'Mortar': {'E': mortar_E, 'v': mortar_v}  # mortar
    }

    # 定义不同相的体积分数
    vol_fractions = {
        'Stone': 0.403,
        'ITZ-2': 0.151,
        'Mortar': 0.423
    }

    # 使用 Mori-Tanaka 方法计算复合材料的有效模量
    K_hom, G_hom, E_hom, v_hom = mori_tanaka_composite(mortar_E, mortar_v, elastic_properties, vol_fractions)
    
    return K_hom, G_hom, E_hom, v_hom

# 处理 pandas 列表中的数据
def calculate_concrete_moduli(df):
    """
    计算 DataFrame 中每一行的体积模量、剪切模量、弹性模量和泊松比，并返回包含计算结果的 DataFrame。
    
    参数：
    - df: 包含 mortar_E 和 mortar_v 的 pandas DataFrame。
    
    返回：
    - 包含 K_bulk_modulus, G_shear_modulus, E_Young_modulus, v_Poisson_ratio 的 DataFrame。
    """
    results = []
    for index, row in df.iterrows():
        mortar_E = row['mortar_E']
        mortar_v = row['mortar_v']
        
        # 调用计算函数
        K_hom, G_hom, E_hom, v_hom = calculate_composite_moduli(mortar_E, mortar_v)
        
        # 保存结果
        results.append([K_hom, G_hom, E_hom, v_hom])
    
    # 返回包含结果的 DataFrame
    return pd.DataFrame(results, columns=['concrete_K', 'concrete_G', 'concrete_E', 'concrete_v'])

# def calculate_compressive_strength(E_values):
#     data = []
    
#     for E in E_values:
#         if pd.isna(E):  # 检查是否为 NaN，如果是则跳过
#             data.append({'Elastic Modulus (GPa)': None, 'Compressive Strength (MPa)': None})
#             continue
        
#         E = float(E)  # 确保 E 是浮点数
#         E_MPa = E * 1000  # 将 GPa 转换为 MPa
#         sigma_B = 10 * (E_MPa / 22000) ** 3  # 使用公式计算抗压强度
#         data.append({'Elastic Modulus (GPa)': E, 'Compressive Strength (MPa)': sigma_B})
    
#     df = pd.DataFrame(data)
    
#     return df[['Compressive Strength (MPa)']]


