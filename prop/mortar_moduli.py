import pandas as pd
import math
import numpy as np

# 将弹性模量和泊松比转换为体积模量和剪切模量



# def E_v_to_K_G(E, v):
#     K = E / (3 * (1 - 2 * v))  # 体积模量 K
#     G = E / (2 * (1 + v))      # 剪切模量 G
#     return K, G

# # 将体积模量和剪切模量转换为弹性模量和泊松比
# def K_G_to_E_v(K, G):
#     E = 9 * K * G / (3 * K + G)  # 弹性模量
#     v = (3 * K - 2 * G) / (2 * (3 * K + G))  # 泊松比
#     return E, v

# # 差分法计算函数
# def differential_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction):
#     """
#     计算基于差分法的复合材料有效模量。
#     """
#     # 定义不同相的弹性参数 (E: 弹性模量, v: 泊松比)
#     elastic_properties = {
#         'Paste': {'E': paste_E, 'v': paste_v},  # 浆体
#         'ITZ-2': {'E': 12.7, 'v': 0.30},       # ITZ-2 (砂浆和砂子之间)
#         'Sand': {'E': 77.6, 'v': 0.15}         # 砂子
#     }

#     # 定义不同相的体积分数
#     vol_fractions = {
#         'Paste': paste_fraction,
#         'ITZ-2': itz_fraction,
#         'Sand': sand_fraction
#     }

#     # 基体 (Paste) 的体积模量和剪切模量
#     K_m, G_m = E_v_to_K_G(paste_E, paste_v)
    
#     # 初始有效模量为基体模量
#     K_eff = K_m
#     G_eff = G_m

#     # 按照差分法逐步引入夹杂物
#     for phase, vol_frac in vol_fractions.items():
#         if phase == 'Paste':  # 跳过基体
#             continue
        
#         # 获取夹杂物的弹性模量和泊松比，并转换为体积模量和剪切模量
#         E_i = elastic_properties[phase]['E']
#         v_i = elastic_properties[phase]['v']
#         K_i, G_i = E_v_to_K_G(E_i, v_i)
        
#         f = vol_frac  # 夹杂物的体积分数
        
#         # 体积模量的差分法公式
#         K_eff = K_eff + f * (K_i - K_eff) / (1 + (K_i - K_eff) / (K_eff + 4 / 3 * G_eff))
        
#         # 剪切模量的差分法公式
#         G_eff = G_eff + f * (G_i - G_eff) / (1 + (G_i - G_eff) / (G_eff + (9 * K_eff + 8 * G_eff) / 6))

#     # 通过体积模量和剪切模量计算弹性模量和泊松比
#     E_eff, v_eff = K_G_to_E_v(K_eff, G_eff)

#     return K_eff, G_eff, E_eff, v_eff

# # 定义函数计算浆体、ITZ和砂子的体积分数   
# def calculate_volume_fractions_with_itz(paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness):
#     """
#     计算浆体、ITZ和砂子的体积分数，基于砂粒的粒径和 ITZ 厚度计算 ITZ 体积。
#     """
#     # 将浆体的ITZ百分比转换为小数
#     itz_fraction_of_paste /= 100
    
#     # 计算浆体的总体积
#     total_paste_volume = paste_weight / paste_density
    
#     # 计算砂子的总体积
#     total_sand_volume = sand_weight / sand_density
    
#     # 计算砂粒的半径和带有 ITZ 的砂粒半径
#     sand_radius = sand_diameter / 2
#     total_radius_with_itz = sand_radius + itz_thickness
    
#     # 计算单个砂粒的体积和带有 ITZ 的砂粒体积
#     single_sand_volume = (4 / 3) * math.pi * (sand_radius ** 3)
#     single_sand_volume_with_itz = (4 / 3) * math.pi * (total_radius_with_itz ** 3)
    
#     # 计算 ITZ 的体积（每个砂粒的 ITZ 体积）
#     itz_volume_per_particle = single_sand_volume_with_itz - single_sand_volume
    
#     # 计算总砂粒数
#     number_of_sand_particles = total_sand_volume / single_sand_volume
    
#     # 计算总的 ITZ 体积
#     total_itz_volume_from_sand = itz_volume_per_particle * number_of_sand_particles
    
#     # 计算浆体中 ITZ 的体积
#     itz_volume_from_paste = total_paste_volume * itz_fraction_of_paste
    
#     # 计算总体系体积（浆体 + 砂子 + ITZ）
#     total_system_volume = total_paste_volume + total_sand_volume + total_itz_volume_from_sand
    
#     # 计算各部分的体积分数
#     paste_fraction = total_paste_volume / total_system_volume
#     itz_fraction = (itz_volume_from_paste + total_itz_volume_from_sand) / total_system_volume
#     sand_fraction = total_sand_volume / total_system_volume
    
#     return paste_fraction, itz_fraction, sand_fraction

# # 计算不同密度、paste_E和paste_v下的弹性模量
# def calculate_mortar_moduli(df):
#     """
#     对pandas DataFrame中的每一行进行计算并返回包含K, G, E, v的DataFrame。
    
#     参数:
#     - df: 包含密度、paste_E和paste_v的pandas DataFrame
#     - paste_weight: 浆体的重量 (g)
#     - itz_fraction_of_paste: ITZ 占浆体体积的百分比
#     - sand_weight: 砂子的重量 (g)
#     - sand_density: 砂子的密度 (g/m³)
#     - sand_diameter: 砂子的粒径 (m)
#     - itz_thickness: ITZ 厚度 (m)
    
#     返回:
#     - 计算结果的pandas DataFrame，包含 K_bulk_modulus, G_shear_modulus, E_Young_modulus, v_Poisson_ratio
#     """
#     paste_weight = 665  # 浆体总重量 (g)
#     itz_fraction_of_paste = 7.3  # ITZ 占浆体体积的百分比
#     sand_weight = 1350  # 砂子的总重量 (g)
#     sand_density = 2570  # 砂子的密度 (g/m³)
#     sand_diameter = 1.5 / 1000  # 砂粒直径 (m)，1.5mm
#     itz_thickness = 20 * 10**-6  # ITZ 厚度 (m)，20微米
#     results = []
#     for index, row in df.iterrows():
#         paste_density = row['density']
#         paste_E = row['E_Young_modulus']
#         paste_v = row['v_Poisson_ratio']
        
#         # 计算体积分数
#         paste_fraction, itz_fraction, sand_fraction = calculate_volume_fractions_with_itz(
#             paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness)
        
#         # 计算四种模量
#         K_hom, G_hom, E_hom, v_hom = differential_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction)
        
#         # 保存结果
#         results.append([K_hom, G_hom, E_hom, v_hom])

#     # 返回包含结果的pandas DataFrame
#     return pd.DataFrame(results, columns=['mortar_K', 'mortar_G', 'mortar_E', 'mortar_v'])

# import numpy as np
# import pandas as pd
# from scipy.optimize import fsolve

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
# def self_consistent_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction):
#     """
#     计算基于自洽方法的复合材料有效模量。
#     """
#     # 定义不同相的弹性参数 (E: 弹性模量, v: 泊松比)
#     elastic_properties = {
#         'Paste': {'E': paste_E, 'v': paste_v},  # 浆体
#         'ITZ-2': {'E': 12.7, 'v': 0.30},       # ITZ-2 (砂浆和砂子之间)
#         'Sand': {'E': 77.6, 'v': 0.15}         # 砂子
#     }

#     # 定义不同相的体积分数
#     vol_fractions = {
#         'Paste': paste_fraction,
#         'ITZ-2': itz_fraction,
#         'Sand': sand_fraction
#     }

#     # 基体 (Paste) 的体积模量和剪切模量
#     K_m, G_m = E_v_to_K_G(paste_E, paste_v)

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

# # 定义函数计算浆体、ITZ和砂子的体积分数   
# def calculate_volume_fractions_with_itz(paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness):
#     """
#     计算浆体、ITZ和砂子的体积分数，基于砂粒的粒径和 ITZ 厚度计算 ITZ 体积。
#     """
#     # 将浆体的ITZ百分比转换为小数
#     itz_fraction_of_paste /= 100
    
#     # 计算浆体的总体积
#     total_paste_volume = paste_weight / paste_density
    
#     # 计算砂子的总体积
#     total_sand_volume = sand_weight / sand_density
    
#     # 计算砂粒的半径和带有 ITZ 的砂粒半径
#     sand_radius = sand_diameter / 2
#     total_radius_with_itz = sand_radius + itz_thickness
    
#     # 计算单个砂粒的体积和带有 ITZ 的砂粒体积
#     single_sand_volume = (4 / 3) * np.pi * (sand_radius ** 3)
#     single_sand_volume_with_itz = (4 / 3) * np.pi * (total_radius_with_itz ** 3)
    
#     # 计算 ITZ 的体积（每个砂粒的 ITZ 体积）
#     itz_volume_per_particle = single_sand_volume_with_itz - single_sand_volume
    
#     # 计算总砂粒数
#     number_of_sand_particles = total_sand_volume / single_sand_volume
    
#     # 计算总的 ITZ 体积
#     total_itz_volume_from_sand = itz_volume_per_particle * number_of_sand_particles
    
#     # 计算浆体中 ITZ 的体积
#     itz_volume_from_paste = total_paste_volume * itz_fraction_of_paste
    
#     # 计算总体系体积（浆体 + 砂子 + ITZ）
#     total_system_volume = total_paste_volume + total_sand_volume + total_itz_volume_from_sand
    
#     # 计算各部分的体积分数
#     paste_fraction = total_paste_volume / total_system_volume
#     itz_fraction = (itz_volume_from_paste + total_itz_volume_from_sand) / total_system_volume
#     sand_fraction = total_sand_volume / total_system_volume
    
#     return paste_fraction, itz_fraction, sand_fraction

# # 计算不同密度、paste_E和paste_v下的弹性模量
# def calculate_mortar_moduli(df):
#     """
#     对pandas DataFrame中的每一行进行计算并返回包含K, G, E, v的DataFrame。
    
#     参数:
#     - df: 包含密度、paste_E和paste_v的pandas DataFrame
    
#     返回:
#     - 计算结果的pandas DataFrame，包含 K_bulk_modulus, G_shear_modulus, E_Young_modulus, v_Poisson_ratio
#     """
#     paste_weight = 665  # 浆体总重量 (g)
#     itz_fraction_of_paste = 7.3  # ITZ 占浆体体积的百分比
#     sand_weight = 1350  # 砂子的总重量 (g)
#     sand_density = 2570  # 砂子的密度 (g/m³)
#     sand_diameter = 1.5 / 1000  # 砂粒直径 (m)，1.5mm
#     itz_thickness = 20 * 10**-6  # ITZ 厚度 (m)，20微米
#     results = []
    
#     for index, row in df.iterrows():
#         paste_density = row['density']
#         paste_E = row['E_Young_modulus']
#         paste_v = row['v_Poisson_ratio']
        
#         # 计算体积分数
#         paste_fraction, itz_fraction, sand_fraction = calculate_volume_fractions_with_itz(
#             paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness)
        
#         # 使用自洽法计算四种模量
#         K_hom, G_hom, E_hom, v_hom = self_consistent_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction)
        
#         # 保存结果
#         results.append([K_hom, G_hom, E_hom, v_hom])

#     # 返回包含结果的pandas DataFrame
#     return pd.DataFrame(results, columns=['mortar_K', 'mortar_G', 'mortar_E', 'mortar_v'])

import numpy as np
import pandas as pd

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
def mori_tanaka_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction):
    """
    计算基于 Mori-Tanaka 方法的复合材料有效模量。
    """
    # 定义不同相的弹性参数 (E: 弹性模量, v: 泊松比)
    elastic_properties = {
        'Paste': {'E': paste_E, 'v': paste_v},  # 基体
        'ITZ-2': {'E': 12.7, 'v': 0.30},       # ITZ-2 (砂浆和砂子之间)
        'Sand': {'E': 77.6, 'v': 0.15}         # 砂子
    }

    # 定义不同相的体积分数
    vol_fractions = {
        'ITZ-2': itz_fraction,
        'Sand': sand_fraction
    }

    # 计算基体 (Paste) 的体积模量和剪切模量
    K_m, G_m = E_v_to_K_G(paste_E, paste_v)

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
        # 获取夹杂物的弹性模量和泊松比，并转换为体积模量和剪切模量
        E_i = elastic_properties[phase]['E']
        v_i = elastic_properties[phase]['v']
        K_i, G_i = E_v_to_K_G(E_i, v_i)

        # 通过 Mori-Tanaka 方法更新基体模量
        K_eff, G_eff = mori_tanaka(K_eff, G_eff, K_i, G_i, vol_frac)

    # 通过体积模量和剪切模量计算弹性模量和泊松比
    E_eff, v_eff = K_G_to_E_v(K_eff, G_eff)

    return K_eff, G_eff, E_eff, v_eff

# 定义函数计算浆体、ITZ和砂子的体积分数   
def calculate_volume_fractions_with_itz(paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness):
    """
    计算浆体、ITZ和砂子的体积分数，基于砂粒的粒径和 ITZ 厚度计算 ITZ 体积。
    """
    # 将浆体的ITZ百分比转换为小数
    itz_fraction_of_paste /= 100
    
    # 计算浆体的总体积
    total_paste_volume = paste_weight / paste_density
    
    # 计算砂子的总体积
    total_sand_volume = sand_weight / sand_density
    
    # 计算砂粒的半径和带有 ITZ 的砂粒半径
    sand_radius = sand_diameter / 2
    total_radius_with_itz = sand_radius + itz_thickness
    
    # 计算单个砂粒的体积和带有 ITZ 的砂粒体积
    single_sand_volume = (4 / 3) * np.pi * (sand_radius ** 3)
    single_sand_volume_with_itz = (4 / 3) * np.pi * (total_radius_with_itz ** 3)
    
    # 计算 ITZ 的体积（每个砂粒的 ITZ 体积）
    itz_volume_per_particle = single_sand_volume_with_itz - single_sand_volume
    
    # 计算总砂粒数
    number_of_sand_particles = total_sand_volume / single_sand_volume
    
    # 计算总的 ITZ 体积
    total_itz_volume_from_sand = itz_volume_per_particle * number_of_sand_particles
    
    # 计算浆体中 ITZ 的体积
    itz_volume_from_paste = total_paste_volume * itz_fraction_of_paste
    
    # 计算总体系体积（浆体 + 砂子 + ITZ）
    total_system_volume = total_paste_volume + total_sand_volume + total_itz_volume_from_sand
    
    # 计算各部分的体积分数
    paste_fraction = total_paste_volume / total_system_volume
    itz_fraction = (itz_volume_from_paste + total_itz_volume_from_sand) / total_system_volume
    sand_fraction = total_sand_volume / total_system_volume
    
    return paste_fraction, itz_fraction, sand_fraction

# 计算不同密度、paste_E和paste_v下的弹性模量
def calculate_mortar_moduli(df):
    """
    对pandas DataFrame中的每一行进行计算并返回包含K, G, E, v的DataFrame。
    """
    paste_weight = 665  # 浆体总重量 (g)
    itz_fraction_of_paste = 7.3  # ITZ 占浆体体积的百分比
    sand_weight = 1350  # 砂子的总重量 (g)
    sand_density = 2570  # 砂子的密度 (g/m³)
    sand_diameter = 1.5 / 1000  # 砂粒直径 (m)，1.5mm
    itz_thickness = 20 * 10**-6  # ITZ 厚度 (m)，20微米
    results = []
    
    for index, row in df.iterrows():
        paste_density = row['density']
        paste_E = row['E_Young_modulus']
        paste_v = row['v_Poisson_ratio']
        
        # 计算体积分数
        paste_fraction, itz_fraction, sand_fraction = calculate_volume_fractions_with_itz(
            paste_weight, paste_density, itz_fraction_of_paste, sand_weight, sand_density, sand_diameter, itz_thickness)
        
        # 使用 Mori-Tanaka 方法计算四种模量
        K_hom, G_hom, E_hom, v_hom = mori_tanaka_composite(paste_E, paste_v, paste_fraction, itz_fraction, sand_fraction)
        
        # 保存结果
        results.append([K_hom, G_hom, E_hom, v_hom])

    # 返回包含结果的pandas DataFrame
    return pd.DataFrame(results, columns=['mortar_K', 'mortar_G', 'mortar_E', 'mortar_v'])

