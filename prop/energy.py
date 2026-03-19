import pandas as pd
import numpy as np

energy_dic={'fly_ash' : 0.033, #in gigajoule per metric ton
            'silica_fume' : 0.036, # values from Alsalaman + al.
            'GGBFS' : 0.857,
            'CSH2' : 1, ##RANDOM VALUE, Please change
            'metakaolin' : 2.5,
            'fine_aggregate' :0.081,
            'coarse_aggregate':0.083,
            'admixture': 29.1,
            'sodium_hydroxide': 20.5,
            'sodium_silicate': 5.371,
            'opc' : 4.53,
            'C3S':4.53,
            'C2S':4.53,
            'C3A':4.53,
            'C4AF':4.53,
            'limestone': 0.02790623,#value from A Review of Carbon Footprint Reduction in Construction Industry, from Design to Operation, Sizirici + al.
            'calcined_clay': 2.734, #best value from Limestone calcined clay cement as a low-carbon solution to meet expanding cement demand in emerging economies, Yudiesky + al.
            'pozzolan': 0.02790623, #value not found (so same as limestone as there are both existing in natural state
           }

emission_dic_old ={#'opc' : 0.74, #CH market ecoinven, 0.84t worldwide average in ton of CO2 (0.73-0.85)
            'C3S':0.2176,
            'C2S':0.1830,
            'C3A':0.1938,
            'C4AF':0.1456,
          'fly_ash' : 0.004,# values from Alsalaman + al
          'CSH2': 0.0082, # RER market ecoinvent
          'silica_fume' : 0.0035, # GLO market ecoinvent 
          'GGBFS' : 0.13, # ROW market ecoinvent
          #'metakaolin' : 0.33, # values from Alsalaman + al
          #'fine_aggregate' :0.0048, # values from Alsalaman + al
          #'coarse_aggregate':0.0048,# values from Alsalaman + al
          #'admixture': 1.88, # values from Alsalaman + al
          #'sodium_hydroxide': 1.915, # values from Alsalaman + al
          #'sodium_silicate': 1.222, # values from Alsalaman + al
          'limestone': 0.0023, # CH market #0.00313,  #value from A Review of Carbon Footprint Reduction in Construction Industry, from Design to Operation  Sizirici + al
          'calcined_clay': 0.27, # ROW market #0.196, #best value from Limestone calcined clay cement as a low-carbon solution to meet expanding cement demand in emerging economies, Yudiesky + al.
          #'pozzolan': 0.00313, #value not found (so same as limestone as there are both existing in natural state
}
emission_dic ={'opc' : 0.74, #CH market ecoinven, 0.84t worldwide average in ton of CO2 (0.73-0.85)
            'C3S':0.82,
            'C2S':0.69,
            'C3A':0.73,
            'C4AF':0.55,
          'fly_ash' : 0.004,# values from Alsalaman + al
          'CSH2': 0.0082, # RER market ecoinvent
          'silica_fume' : 0.0035, # GLO market ecoinvent 
          'GGBFS' : 0.13, # ROW market ecoinvent
          'metakaolin' : 0.33, # values from Alsalaman + al
          'fine_aggregate' :0.0048, # values from Alsalaman + al
          'coarse_aggregate':0.0048,# values from Alsalaman + al
          'admixture': 1.88, # values from Alsalaman + al
          'sodium_hydroxide': 1.915, # values from Alsalaman + al
          'sodium_silicate': 1.222, # values from Alsalaman + al
          'limestone': 0.0023, # CH market #0.00313,  #value from A Review of Carbon Footprint Reduction in Construction Industry, from Design to Operation  Sizirici + al
          'calcined_clay': 0.27, # ROW market #0.196, #best value from Limestone calcined clay cement as a low-carbon solution to meet expanding cement demand in emerging economies, Yudiesky + al.
          'pozzolan': 0.00313, #value not found (so same as limestone as there are both existing in natural state
}


# 

# def energy_emission_data(data):
#     ''' Compute the energy need and the Co2 emission for a given recipe'''

#     # 过滤掉不需要的成分 (例如 'opc')
#     filtered_data = data.loc[:, energy_dic.keys() & emission_dic.keys()]

#     # 计算每列的总能耗
#     energy_needs = (filtered_data * pd.Series(energy_dic)).sum(axis=1)
    
#     # 计算每列的总排放
#     emissions = (filtered_data * pd.Series(emission_dic)).sum(axis=1)

#     # 计算每列的总材料重量
#     material_weights = filtered_data.sum(axis=1)

#     # 避免除以零的情况
#     with np.errstate(divide='ignore', invalid='ignore'):
#         energy_per_kg = energy_needs / material_weights
#         energy_per_kg[material_weights == 0] = np.nan  # 处理除以零的情况

#         emission_per_kg = emissions / material_weights
#         emission_per_kg[material_weights == 0] = np.nan  # 处理除以零的情况

#     # 创建最终的 DataFrame
#     energy_emission_results = pd.DataFrame({
#         'Energy_per_kg': energy_per_kg,
#         'Emission_per_kg': emission_per_kg
#     })

#     return energy_emission_results

# def energy_emission_data(data):
#     ''' Compute the energy need and the Co2 emission for a given recipe'''
#     energy_np = np.full(data.shape[1], 100*energy_dic['opc']) 
#     emission_np = np.full(data.shape[1], 100*emission_dic['opc'])
#     mat_quant_np = np.full(data.shape[1], 100)
#     for i in range(data.shape[1]):
#         dset = data.iloc[:,i]
#         for key in emission_dic.keys():
#             if key in data.index:
#                 energy_np[i] += dset.loc[key]*energy_dic[key]
#                 emission_np[i] += dset.loc[key]*emission_dic[key]
#                 mat_quant_np[i] += dset.loc[key]
#     return pd.DataFrame(energy_np/mat_quant_np), pd.DataFrame(emission_np/mat_quant_np)


def CO2_emission(data):
    
    # Filter data to only those ingredients that have emission factors
    filtered_data = data.loc[:,emission_dic.keys()]

    # Compute total CO2 emissions per column
    co2_emissions = (filtered_data * pd.Series(emission_dic)).sum(axis = 1)

    # Compute total material weights per column
    material_weights = filtered_data.sum(axis = 1)

    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        co2_per_kg = co2_emissions / material_weights
        co2_per_kg[material_weights == 0] = np.nan  # Handle division by zero if any

    # Creating the final DataFrame
    total_CO2_emissions = pd.DataFrame({
        #'CO2': co2_emissions,
        #'weight': material_weights,
        'CO2_per_kg': co2_per_kg
    })

    return total_CO2_emissions

# def energy_emission_data(data):
#     # Filter data to only those ingredients that have emission factors
#     filtered_data = data.loc[:,energy_dic.keys()]

#     # Compute total CO2 emissions per column
#     cement_energy = (filtered_data * pd.Series(energy_dic)).sum(axis = 1)

#     # Compute total material weights per column
#     material_weights = filtered_data.sum(axis = 1)

#     # Avoid division by zero
#     with np.errstate(divide='ignore', invalid='ignore'):
#         cement_per_kg = cement_energy / material_weights
#         cement_per_kg[material_weights == 0] = np.nan  # Handle division by zero if any

#     # Creating the final DataFrame
#     total_cement_energy = pd.DataFrame({
#         #'CO2': co2_emissions,
#         #'weight': material_weights,
#         'cement_per_kg': cement_per_kg
#     })

#     return total_cement_energy
def energy_emission_data(data):
    ''' Compute the energy need and the Co2 emission for a given recipe'''

    # 找出既在 data 列中又在 energy_dic 中的成分
    common_ingredients = list(set(data.columns).intersection(set(energy_dic.keys())))
    
    # 过滤掉不存在的列
    filtered_data = data.loc[:, common_ingredients]

    # 计算每列的总能耗
    energy_needs = (filtered_data * pd.Series(energy_dic)).sum(axis=1)
    
    # 计算每列的总排放
    emissions = (filtered_data * pd.Series(emission_dic)).sum(axis=1)

    # 计算每列的总材料重量
    material_weights = filtered_data.sum(axis=1)

    # 避免除以零的情况
    with np.errstate(divide='ignore', invalid='ignore'):
        energy_per_kg = energy_needs / material_weights
        energy_per_kg[material_weights == 0] = np.nan  # 处理除以零的情况

        emission_per_kg = emissions / material_weights
        emission_per_kg[material_weights == 0] = np.nan  # 处理除以零的情况

    # 创建最终的 DataFrame
    energy_emission_results = pd.DataFrame({
        'Energy_per_kg': energy_per_kg,
        'Emission_per_kg': emission_per_kg
    })

    return energy_emission_results

def emission(x):
    ''' Compute the energy need and the Co2 emission for a given recipe'''
    emission_np = np.full(data.shape[1], 100*emission_dic['opc'])
    mat_quant_np = np.full(data.shape[1], 100)
    for i in range(data.shape[1]):
        dset = data.iloc[:,i]
        for key in emission_dic.keys():
            if key in data.index:
                energy_np[i] += dset.loc[key]*energy_dic[key]
                print (dset.loc[key]*energy_dic[key])
                emission_np[i] += dset.loc[key]*emission_dic[key]
                mat_quant_np[i] += dset.loc[key]
    return pd.DataFrame(energy_np/mat_quant_np), pd.DataFrame(emission_np/mat_quant_np)