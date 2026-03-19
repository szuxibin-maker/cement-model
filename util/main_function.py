import numpy as np
from run.hydration import run_hydration
from run.hydration import parrot_killoh
from run.hydration import to_phase_first_dict
from run.hydration import plot_bars
from run.hydration import phase_plot
from run.hydration import co2_phase_plot
from util.final_hydration import hydration

def plot_result(data_recipe, output_times):
    # 初始化存储结果的列表
    pH_dic_list = []
    density_final_list = []

    # 设置常量
    recipe_meas = data_recipe

    # 进行水化过程并获取结果
    index_to_drop, gems_vol, gems_phase_masses, density, gems_pH_dic, gems_B, element_CSHQ, gems_CSH_species_mass, gems_CSH_species_volumes, element_aq = hydration(recipe_meas, output_times, fail=False)
    
    if gems_vol is None:
        print("GEMS volume is None, corresponding recipe_meas:", recipe_meas)
        return index_to_drop, None, None, None, None, None, None, None, None
    # print(index_to_drop)
    classfication_CSH = classify_CSH_species(gems_CSH_species_volumes)
    update_gems_vol_with_CSH(gems_vol, classfication_CSH)

    # 处理体积分数
    gems_vol_phase_first_inital = to_phase_first_dict(gems_vol)
    # 计算基于0时刻总体积的体积分数
    gems_vol_frac_based_on_time_0 = calculate_volume_fraction(gems_vol_phase_first_inital)
    
    # 计算基于各自时间点总体积的体积分数
    gems_vol_frac_based_on_each_time = calculate_volume_fraction_based_on_each_time(gems_vol_phase_first_inital, output_times)
    
    # 处理质量数据
    classfication_CSH_mass = classify_CSH_mass(gems_CSH_species_mass)
    update_gems_phase_masses_with_CSH(gems_phase_masses, classfication_CSH_mass)
    gems_masses_phase_first = to_phase_first_dict(gems_phase_masses)

    # 创建 gems_vol_phase_first_inital 的副本
    gems_vol_phase_first_inital_delet_gas = gems_vol_phase_first_inital.copy()
    
    # 删除 'gas_gen' 键
    del gems_vol_phase_first_inital_delet_gas['gas_gen']
    
    # 后续计算
    total_mass_initial = sum(values[0] for values in gems_masses_phase_first.values())
    result_1 = sum_values_at_positions(gems_vol_phase_first_inital_delet_gas)
    density_final = [(total_mass_initial * 0.001) / vol for vol in result_1]
    density_final = {output_times[i]: density_final[i] for i in range(len(output_times))}
    
    # 存储结果
    density_final_list.append(density_final)
    pH_dic_list.append(gems_pH_dic)
    print('gems_B',gems_B)
    # 返回计算结果
    return index_to_drop, gems_vol_frac_based_on_each_time, gems_vol_frac_based_on_time_0, gems_masses_phase_first, density_final_list, pH_dic_list, gems_B, element_CSHQ, element_aq

def classify_CSH_species(gems_CSH_species_volumes):
    """根据体积分类CSH物质"""
    classfication_CSH = {}
    for time, data in gems_CSH_species_volumes.items():
        JenD=data['CSHQ-JenD']
        JenH=data['CSHQ-JenH']
        TobD=data['CSHQ-TobD']
        TobH = data['CSHQ-TobH'] #+ data['KSiOH'] + data['NaSiOH']
        classfication_CSH[time] = {'CSHQ-JenD': JenD,
              'CSHQ-JenH':JenH,
               'CSHQ-TobD':TobD,                    
            'CSHQ-TobH':TobH,
        }
    return classfication_CSH

def classify_CSH_mass(gems_CSH_species_mass):
    """根据质量分类CSH物质"""
    classfication_CSH_mass = {}
    for time, data in gems_CSH_species_mass.items():
        JenD=data['CSHQ-JenD']
        JenH=data['CSHQ-JenH']
        TobD=data['CSHQ-TobD']
        TobH = data['CSHQ-TobH'] #+ data['KSiOH'] + data['NaSiOH']
        classfication_CSH_mass[time] = {'CSHQ-JenD':JenD,
              'CSHQ-JenH':JenH,
               'CSHQ-TobD':TobD,                    
            'CSHQ-TobH':TobH,
        }
    return classfication_CSH_mass

def update_gems_vol_with_CSH(gems_vol, classfication_CSH):
    """更新gems_vol字典，将CSH分类信息合并进去"""
    for time in gems_vol:
        if time in classfication_CSH:
            gems_vol[time].update(classfication_CSH[time])
        if 'CSHQ' in gems_vol[time]:
            del gems_vol[time]['CSHQ']

def update_gems_phase_masses_with_CSH(gems_phase_masses, classfication_CSH_mass):
    """更新gems_phase_masses字典，将CSH质量分类信息合并进去"""
    for time in gems_phase_masses:
        if time in classfication_CSH_mass:
            gems_phase_masses[time].update(classfication_CSH_mass[time])
        if 'CSHQ' in gems_phase_masses[time]:
            del gems_phase_masses[time]['CSHQ']

def calculate_volume_fraction(gems_vol_phase_first_inital):
    """计算体积分数，排除掉气体相和水相"""
    gems_vol_frac_phase_first_gas_gen = {key: value for key, value in gems_vol_phase_first_inital.items() if key not in ['gas_gen']}
    total_delet_gas_gen = sum(values[0] for values in gems_vol_frac_phase_first_gas_gen.values())
    gems_vol_frac_phase_first = {key: values / total_delet_gas_gen for key, values in gems_vol_frac_phase_first_gas_gen.items()}
    return gems_vol_frac_phase_first

def calculate_volume_fraction_based_on_each_time(gems_vol_phase_first_inital, output_times):
    """计算基于每个时间点总体积的体积分数，剔除掉气体相，并返回一个包含所有时间点的 NumPy 数组"""
    
    # 初始化一个空的字典来存储体积分数
    gems_vol_frac_based_on_each_time = {key: np.zeros(len(output_times)) for key in gems_vol_phase_first_inital.keys() if key != 'gas_gen'}
    
    for time_index, time in enumerate(output_times):
        # 获取当前时间点的体积，并剔除 gas_gen 相
        phase_volumes = {key: values[time_index] for key, values in gems_vol_phase_first_inital.items() if key != 'gas_gen'}
        
        # 计算该时间点剔除 gas_gen 后的总和
        total_volume_excluding_gas = sum(phase_volumes.values())
        
        # 计算该时间点下每个相的体积分数并存储在数组中
        for key in gems_vol_frac_based_on_each_time.keys():
            gems_vol_frac_based_on_each_time[key][time_index] = phase_volumes[key] / total_volume_excluding_gas
    
    return gems_vol_frac_based_on_each_time


def porosity(vol_frac_phase):
   # 删除 'aq_gen' 项目
    filtered_data = {key: value for key, value in vol_frac_phase.items() if key != 'aq_gen'}

# 计算所有剩余项目的逐元素相加
    sum_all_phase = np.sum(np.array(list(filtered_data.values())), axis=0)

# 计算 1 - sum_all_phase 得到 result
    result = 1 - sum_all_phase
    return result

def sum_values_at_positions(data_dict):
    # 获取所有的值（数组）
    arrays = list(data_dict.values())
    sum_array = np.sum(np.array(arrays), axis=0)
    
    return sum_array

