import pandas as pd
import numpy as np


CaO_in_phase = {'CSHQ-JenD': 0.497, 'CSHQ-JenH': 0.43 ,'CSHQ-TobD': 0.39,'CSHQ-TobH': 0.3,
                #TobD , JenH JenD s三者平均值(, # K  obtain in Modeling the elastic properties of Portland cement paste Haecker+al.有两种一种是高，一种是低
            'C3(AF)S0_84H':0.323,#Ca3Fe2(SiO4)0.84 (OH)8.64
            'C3FS1_34H3_32':0.3593,#Ca3Fe2(SiO4)1.34 (OH)6.64
            'C4Ac0_5H12':0.322,#Ca4Al2(CO3)0.5(OH)13·7H2O
            'C4AsH14':0.312,#Ca4Al2(SO4)(OH)12·6H2O
            'CO3_SO4_AFt':0.27,#tricarboalu03:(CO3)Ca2Al0.6666667(OH)4(H2O)8.6666667,M:382.32,0.2934 , ettringite03_ss:(SO4)Ca2Al0.6666667(OH)4(H2O)8.6666667,M:418.3703 :0.2681
            'SO4_CO3_AFt':0.27,##tricarboalu03:(CO3)Ca2Al0.6666667(OH)4(H2O)8.6666667,M:382.32,0.2934 , ettringite03_ss:(SO4)Ca2Al0.6666667(OH)4(H2O)8.6666667,M:418.3703 :0.2681
            'OH_SO4_AFm':0.34,#C4AH13:,Ca4Al2(OH)14(H2O)6,M:560.5,0.334,  monosulphate12:Ca4Al2SO10(H2O)12,M:622.52,0.3603
            'Portlandite' : 0.757,#
            'ettringite':0.27,#ettringite: ((H2O)2)Ca6Al2(SO4)3(OH)12(H2O)24,M:1255.11, 0.2681,ettringite30 ,Ca6Al2(SO4)3(OH)12(H2O)24,M:1219.08, 0.276
            'ettringite-FeAl':0.2681,#ettringite05:Ca3Al(SO4)1.5(OH)6(H2O)13,M:627.5554,0.2681, Fe-ettringite05 :,Ca3Fe|3|(SO4)1.5(OH)6(H2O)13,M:656.4,0.219
             'ettringite-AlFe':0.2681,#ettringite05:Ca3Al(SO4)1.5(OH)6(H2O)13,M:627.5554,0.2681, Fe-ettringite05 :,Ca3Fe|3|(SO4)1.5(OH)6(H2O)13,M:656.4,0.219
            'straetlingite': 0.27,#straetlingite:Ca2Al2SiO7(H2O)8,M:418.3228,0.2681,straetlingite7:Ca2Al2SiO7(H2O)7,M:400.3,0.2802
            'C4AcH11':0.3946}#4CaO⋅Al 

# 'C3(AF)S0.84H','C3FS1.34H3.32','C4Ac0.5H12', 'C4AsH14', 'CO3_SO4_AFt','OH-hydrotalcite', 'OH_SO4_AFm', 'Portlandite','ettringite', 'ettringite-FeAl','straetlingite'
# ,
               # 'HCA-Friedels-ss': 0.4 }#AFt}'HCA-Friedels-ss': 0.4

# 创建质量的字典，不进行归一化
def create_mass_phase(index, data):
    mass_phase_dic = {index[i]: data.iloc[i] for i in range(index.shape[0])}
    return mass_phase_dic

# 计算CO2吸收量
def compute_CO2_absorption(mass_phase):
    ''' 计算CO2 吸收量 '''
    total_CaO = 0
    for key, value in mass_phase.items():
        if key in CaO_in_phase:
            # 质量乘以相中的CaO含量
            total_CaO += value * CaO_in_phase[key]
    # CaO 总量乘以0.785计算CO2吸收量
    CO2_absorption = total_CaO * 0.785
    return CO2_absorption

# 数据处理和CO2吸收量计算
def data_CO2_absorption(data):
    ''' 计算数据集中的 CO2 吸收量 '''
    global mass_phase
    CO2_absorption_np = np.zeros(data.shape[1])
    index = data.index
    for i in range(data.shape[1]):
        mass_phase = create_mass_phase(index, data.iloc[:, i])
        CO2_absorption = compute_CO2_absorption(mass_phase)
        CO2_absorption_np[i] = CO2_absorption
    return pd.DataFrame(CO2_absorption_np, columns=['CO2_absorption'])



import pandas as pd

# 定义每种材料中的 CaO 含量字典
CaO_in_reactant_dic = {
    'fly_ash': 0.0991,
    'silica_fume': 0.0015,
    'GGBFS': 0.41,
    'CSH2': 0.298,
    'metakaolin': 0,
    'fine_aggregate': 0,
    'coarse_aggregate': 0,
    'admixture': 0,
    'sodium_hydroxide': 0,
    'sodium_silicate': 0,
    'opc': 0,
    'C3S': 0.737,
    'C2S': 0.651,
    'C3A': 0.623,
    'C4AF': 0.461,
    'limestone': 0.423,
    'calcined_clay': 0.52,
    'pozzolan': 0.0905
}

# 计算每一行的总 CaO 含量
def calculate_total_CaO_in_reactant(data: pd.DataFrame) -> pd.DataFrame:
    """
    根据数据集中每一行的成分数量计算总 CaO 含量。

    参数:
    - data (pd.DataFrame): 数据集，其中包含各成分的数量，列名应与 CaO_in_reactant_dic 的键匹配。

    返回:
    - total_CaO_df (pd.DataFrame): 包含每行总 CaO 含量的 DataFrame。
    """
    # 保留数据集中同时存在于 CaO_in_reactant_dic 中的成分
    common_columns = list(set(data.columns).intersection(set(CaO_in_reactant_dic.keys())))
    filtered_data = data[common_columns]
    
    # 计算每一行的 CaO 含量
    total_CaO_per_row = (filtered_data * pd.Series(CaO_in_reactant_dic)).sum(axis=1)
    
    # 返回结果作为 DataFrame
    return pd.DataFrame(total_CaO_per_row, columns=['Total_CaO_in_reactant'])

import pandas as pd

# 定义每种材料中的 CaO 含量字典
CaO_in_reactant_dic_time = {
    'FA': 0.0991,
    'SF': 0.0015,
    'Slag': 0.41,
   'Gp': 0.298,
    'metakaolin': 0,
    'fine_aggregate': 0,
    'coarse_aggregate': 0,
    'admixture': 0,
    'sodium_hydroxide': 0,
    'sodium_silicate': 0,
    'opc': 0,
    'Alite': 0.737,
    'Belite': 0.651,
    'Aluminate': 0.623,
    'Ferrite': 0.461,
    'limestone': 0.423,
    'clay': 0.52,
    'pozzolan': 0.0905
}

# 计算每一行的总 CaO 含量
def calculate_total_CaO_in_reactant_time(data: pd.DataFrame) -> pd.DataFrame:
    """
    根据数据集中每一行的成分数量计算总 CaO 含量。

    参数:
    - data (pd.DataFrame): 数据集，其中包含各成分的数量，列名应与 CaO_in_reactant_dic 的键匹配。

    返回:
    - total_CaO_df (pd.DataFrame): 包含每行总 CaO 含量的 DataFrame。
    """
    # 保留数据集中同时存在于 CaO_in_reactant_dic 中的成分
    common_columns = list(set(data.columns).intersection(set(CaO_in_reactant_dic_time.keys())))
    filtered_data = data[common_columns]
    
    # 计算每一行的 CaO 含量
    total_CaO_per_row = (filtered_data * pd.Series(CaO_in_reactant_dic_time)).sum(axis=1)
    
    # 返回结果作为 DataFrame
    return pd.DataFrame(total_CaO_per_row, columns=['Total_CaO_in_reactant_time'])
