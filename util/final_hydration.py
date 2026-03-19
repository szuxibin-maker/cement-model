from run.GEMSCalc import GEMS
from run.hydration import parrot_killoh
import pandas as pd
import numpy as np
import cmocean

# dictionary for the recipe of some SCM from CEMGEMS CEM-II-BV based on 100 g of each SCM  always give formula in moles!!!
SCM = { 'limestone': ({'Al' : 0.052961280988, 'C' : 0.85663128045, 'Ca': 0.75431457236,'Fe' : 0.025048813876, 'K': 0.012739394454,'Mg' : 0.044660136362, 'Na': 0.016134497168, 'O': 3.0564427725, 'Si':0.20637670739},  'kg'),
        'silica_fume':({'Al' : 0.005883425509, 'Ca': 0.002674346509,'Fe' : 0.0002504387073, 'K': 0.0063684402277,'Mg' :0.00099225162461, 'Mn': 0.0001266581442,'Na' : 0.0029036363721, 'O': 3.3128876396,'P' : 0.00028174463078, 'Si': 1.6423635207,'Ti' : 0.0048814082984}, 'kg'),
      'GGBFS': ({'Al' : 0.25068590354, 'C' : 0.026137740473, 'Ca': 0.7400381436, 'Fe': 0.0062268216656, 'H': 0.24594931045, 'K': 0.0066234897348, 'Mg':0.18043867977, 'Na': 0.0050333680703, 'O':2.568871619, 'S': 0.034432524341, 'Si': 0.54806878508, 'Ti':0.0064682389664},'kg'), 
       'calcined_clay':({'Al' : 0.86001188208, 'Ca' : 0.0017850310446, 'Fe': 0.0037610756424, 'H': 0.16669147874, 'K': 0.0021253535161, 'Na':0.009690369309, 'O': 3.1677258055, 'P':0.0028208175955, 'S':0.0012502298115, 'Si':0.86631529281, 'Ti': 0.018797209003},'kg'), 
       'fly_ash': ({'Al' : 0.35503673551, 'C' : 0.016132843743, 'Ca': 0.17672003338, 'Fe': 0.11522454383, 'K': 0.051594547539, 'Mg':0.087831601512, 'Na': 0.03259168428, 'O':2.8800652865,'S': 0.0049959283184, 'Si':  0.91038757213},'kg'),
      'Pozzolan': ({'Al' : 0.3500773, 'Ca' : 0.1609977,  'Fe': 0.05360106, 'H':0.3377904, 'K': 0.1611911, 'Mg':0.03044469, 'Na': 0.09915054, 'O':2.902273,'S': 0.008098946, 'Si':  0.8813103},'kg'),
       'biochar': ({'Al' : 0.01473923, 'Ca' : 1.076925,  'Fe': 0.03687062,  'K': 0.4672641, 'Mg':0.1200335, 'Na': 0.08669282, 'O':1.781338,'P':0.02683121,'S': 0.04139704, 'Si':  0.01935863},'kg'),
        # 'silica_fume': ({'Al' :  0.070, 'C' : 0.03, 'Ca' :  1.2 ,  'Fe':  0.003 ,  'K': 0.009, 'Mg': 0.027, 'Na': 0.001, 'O':2.228,'P':0.0005,'S': 0.042 , 'Si':  0.363,'Ti' : 0.0001,'P' :0.001},'mol')#这个是wpc
      }
SCM_keys = ['limestone', 'silica_fume', 'GGBFS', 'calcined_clay', 'fly_ash', 'Pozzolan', 'biochar']


'''
Formula in %mass for 100g: 
- Limestone: CaO: 42.3, Al2O3: 2.7, Fe2O3: 2, SiO2: 12.4, Na2O: 0.5, MgO: 1.8, K2O: 0.6, CO2: 37.7

- silica_fume: CaO: 0.14997, MgO: 0.039992, K2O: 0.29994, Na2O: 0.089982, TiO2: 0.38992, Mn2O3: 0.009998, P2O5: 0.019996, Al2O3: 0.29994, Fe2O3: 0.019996, SiO2: 98.68

- GGBFS: SiO2: 32.93, Al2O3: 12.78, Fe2O3: 0.49717, CaO: 41.499, MgO: 7.2724, SO3: 0.14623, Na2O: 0.15598, K2O: 0.31195, H2S: 1.1113, H2O: 1.628, TiO2: 0.51667, CO2: 1.1503

- calcined clay: SiO2: 52.052, Al2O3: 43.844, Fe2O3: 0.3003, CaO: 0.1001, SO3: 0.1001, Na2O: 0.3003, K2O: 0.1001, TiO2: 1.5015, P2O5: 0.2002, H2O: 1.5015 

- fly ash: Ca0: 9.91, Al2O3: 18.1, Fe2O2: 9.2, SiO2: 54.7, CO2: 0.71, MgO: 3.54, K2O: 2.43, Na2O: 1.01, SO2: 0.4, P2O5: 0, C: 0
-Pozzolan ：CaO: 9.05, Al2O3: 17.89, Fe2O3: 4.29, SiO2: 53.08, H20:3.05 , MgO: 1.23, K2O: 7.61, Na2O: 3.08, SO3: 0.65, TiO2: 0.31
- biochar (forest wood): CaO: 58.67, Al2O3: 0.73, Fe2O3: 2.86, SiO2: 1.13, MgO: 4.7, K2O: 21.38, Na2O: 2.61, SO3: 3.22,Cl:0.47,P2O5: 1.85
 reference: Effect of cement partial substitution by waste-based biochar in mortars properties
'''
'''
ettringite:Ca6Al2(SO4)3(OH)12(H2O)26
monosulfate: Ca4Al2SO10(H2O)12
hydrotalcite:Mg4Al2O7(H2O)10

'''



element_names = [
"Al", "B", "C", "Ca", "Cl", "Fe", "H", "K", "Mg", "Mn", 
        "Na", "Nit", "O", "P", "S", "Si", "Ti", "Zz", "OH-"]



# dictionary for the 5PL model (values found in cemgem)
dic_5lp = {'GGBFS' : {'A':0, 'B':1, 'C':7, 'D':70, 'G':1},#{'A':0, 'B':1, 'C':7, 'D':70, 'G':1},
           'calcined_clay' : {'A':0, 'B':1, 'C':5, 'D':50, 'G':1},
           'fly_ash' : {'A':0, 'B':0.7, 'C':85.1, 'D':60, 'G':1},
          'Pozzolan' : {'A':0, 'B':0.7, 'C':85.1, 'D':0, 'G':0},
          'silica_fume' : {'A':0.014, 'B':1, 'C':16.62, 'D':74.98, 'G':1},# (34 + 10 * time)−0.26 REVERSE FROM Hydration of a silica fume blended low-alkali shotcrete cement
          'limestone' : {'A':100 ,'B':1, 'C':50, 'D':100, 'G':1},
          'biochar': {'A':0.00014, 'B':1, 'C':16.62, 'D':0.7498, 'G':1}}#FULL REACTION DEGREE,Pozzolan,silica_fume,limestone

# def silica_fume_hydration ():
    

# 各种材料的5LP参数（来自cemgem）
DIC_5LP = {
    'GGBFS': {'A': 0, 'B': 1, 'C': 7, 'D': 70, 'G': 1},
    'calcined_clay': {'A': 0, 'B': 1, 'C': 5, 'D': 50, 'G': 1},
    'fly_ash': {'A': 0, 'B': 0.7, 'C': 85.1, 'D': 60, 'G': 1},
    'Pozzolan': {'A': 0, 'B': 0.7, 'C': 85.1, 'D': 0, 'G': 0},
    'silica_fume': {'A': 0.014, 'B': 1, 'C': 16.62, 'D': 74.98, 'G': 1},
    'limestone': {'A': 100, 'B': 1, 'C': 50, 'D': 100, 'G': 1},
    'OPC': {'A': 2.06, 'B': 1.38, 'C': 4.98, 'D': 92.4, 'G': 0.5},
    # 别名映射
    # 'Slag': {'A': 0, 'B': 1, 'C': 7, 'D': 70, 'G': 1},  # same as GGBFS
    # 'Fly Ash': {'A': 0, 'B': 0.7, 'C': 85.1, 'D': 60, 'G': 1},  # same as fly_ash
    # 'Silica Fume': {'A': 0.014, 'B': 1, 'C': 16.62, 'D': 74.98, 'G': 1},  # same as silica_fume
    # 'PC': {'A': 2.06, 'B': 1.38, 'C': 4.98, 'D': 92.4, 'G': 0.5},  # same as OPC
}
SCM_density = {
    'GGBFS': 2880,
    'calcined_clay': 2300,
    'fly_ash': 2490,
    'silica_fume': 3200,
    'Pozzolan': 2000,
    'limestone': 2710,
    'Pozzolan': 2000,
}
clinker_density = {
    'C3S': 3120,
    'C2S': 3326,
    'C3A': 3029,
    'C4AF': 3732,
    'CSH2': 2305,
}

def add_material_to_gemsk(gemsk, recipe):
    ''' Add materials from recipe to gemsk object '''

    all_species = recipe['clink_phases'].copy()
    all_species["H2O@"] = recipe['wc']
    all_species["CO2"] = recipe['CO2']
    all_species["Gp"]=recipe['CSH2']
    # print(all_species)
    # print('recipe',recipe)
    for name in all_species.keys():
        all_species[name]*=1e-3
    gemsk.T = recipe['T'] + 273.15
    gemsk.add_multiple_species_amt(all_species,units = "kg")
    for component in recipe.keys(): # add SCM
        if component not in ['clink_phases', 'wc', 'CSH2', 'T', 'RH', 'fineness','CO2']:
            #print(component)
            assert component in SCM.keys()
            value = recipe[component] *1e-3   # we give the input recipe in g
            gemsk.add_amt_from_formula(SCM[component][0], value, units='kg')  
            # print(SCM[component][0])


    gemsk.add_species_amt("O2",1e-6) # to reduce stiffness related to redox

def add_material_for_hydration_to_gemsk(gemsk, recipe):
    ''' Add materials from recipe to gemsk object '''
    all_species = recipe['clink_phases'].copy()
    all_species["H2O@"] = recipe['wc'] 
    all_species["CO2"] = recipe['CO2']
    all_species["Gp"]=recipe['CSH2']
    for name in all_species.keys():
        all_species[name]*=1e-3
    gemsk.T = recipe['T'] + 273.15
    gemsk.add_multiple_species_amt(all_species,units = "kg")
    gemsk.add_species_amt("O2",1e-6) # to reduce stiffness related to redox

def process_gems_result(gemsk):
    ''' Function to process gem result once the equilibration is done '''
    # gemsk.supress_phase('gas_gen')
    gems_vol=gemsk.phase_volumes
    # print('gems_vol',gems_vol)
    total_vol_delet_gas_sum = sum(value for key, value in gems_vol.items() if key != 'gas_gen')
    other_vol_delete_gas={key: value for key, value in gems_vol.items() if key not in [ 'gas_gen']}
    # print('other_vol_delete_gas',other_vol_delete_gas)
    vol_delete_gas=other_vol_delete_gas.copy()
    for key,val in vol_delete_gas.items():
        vol_delete_gas[key]/=total_vol_delet_gas_sum
    
    gems_vol_frac = vol_delete_gas
    # print('gems_vol_frac',gems_vol_frac)
    gems_phase_amounts = gemsk.phase_masses #.phase_masses是质量，phase_amounts是mol
    gems_B = gemsk.b
    gems_pE = gemsk.pE
    gems_pH = gemsk.pH
    element_CSHQ=gemsk.reset_CSHQ_composition
    element_aq=gemsk.aq_element
  
    element_CSHQ = dict(zip(element_names, element_CSHQ))
    element_aq = dict(zip(element_names, element_aq))
    total_mass_delet_gas_gen=sum(value for key, value in (gemsk.phase_masses).items() if key != 'gas_gen')
    gems_density = total_mass_delet_gas_gen/total_vol_delet_gas_sum
    print('total_mass_delet_gas_gen',total_mass_delet_gas_gen)
    gems_aq_volume_frac=gemsk.aq_volume_frac
    # print('other_vol_delete_gas',other_vol_delete_gas)
    other_vol_delete_gas = {key: value / 0.000001 for key, value in other_vol_delete_gas.items()}
    return gems_vol_frac , gems_phase_amounts, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac,element_CSHQ,element_aq,other_vol_delete_gas



    
def final_hydration(recipe, fail = False):
    ''' Function to compute the final hydration from a recipe '''
    input_file = 'gems_files/MySystem2-dat.lst'
    gemsk = GEMS(input_file)
    # del recipe['wc_ratio']
    add_material_to_gemsk(gemsk, recipe)
     # === AFt-phases（三钙铝石类，AFt） ===
    # gemsk.supress_phase("ettringite-AlFe")     # 含铝铁的钙矾石，AFt 结构变体之一 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    # gemsk.supress_phase("ettringite-FeAl")     # 含铁铝的钙矾石，另一种AFt结构混合型 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    # gemsk.supress_phase("ettringite")          # 钙矾石，C6AŜ3H32，典型的硫铝酸盐水化产物
    # # gemsk.supress_phase("SO4_CO3_AFt")         # 硫酸根-碳酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    # gemsk.supress_phase("CO3_SO4_AFt")         # 碳酸根-硫酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    # gemsk.supress_phase("C6AsH13")             # 硫酸型AFt矿物变体，可能水合度不同 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    # gemsk.supress_phase("C6AsH9")              # 更低水合度的硫酸型AFt矿物 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    gemsk.supress_phase("thaumasite")          # 硅钙矾石（Thaumasite），Si 替代 SO₄ 形成的 AFt 型结构 需低温（<15°C）+高 CO₂ + Si源（如硅灰）+ SO₄²⁻，在 OPC 正常条件下极少生成，但在寒冷潮湿环境下可能导致腐蚀问题（称为“thaumasite 硫酸盐攻击”）
    
    # === AFm-Silica（硅铝酸盐型 AFm） ===
    # gemsk.supress_phase("straetlingite")       # 含硅 AFm，相当于 C2ASH8，常见于硅含量较高体系 在硅铝材料体系（如含火山灰、偏高岭土的胶凝材料）中可见，在纯 OPC 中较少见，但并非罕见
    
    # === AFm-phases（单盐型水化铝酸盐 AFm） ===
    # gemsk.supress_phase("C4AH19")              # 高水合度铝酸钙水合物，早期水化产物 高水合度，早期生成，
    # gemsk.supress_phase("C4AH13")              # 中等水合度 AFm，常在水化过程中生成 ，是常见稳定产物之一
    # gemsk.supress_phase("C4AH11")              # 较低水合度的 AFm，相对更稳定 ，相对稳定，晚期可能形成
    # gemsk.supress_phase("CAH10")               # 一种低温条件下形成的 AFm 类型 低温环境或高水灰比条件下更易生成，实验中难观察
    # gemsk.supress_phase("C4Ac0.5H12")          # 碳酸-氢氧共插层型 AFm，相对较高水合度（CO₃²⁻ + OH⁻）需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4Ac0.5H105")         # 同上，略低水合度版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4Ac0.5H9")           # 碳酸型 AFm 的低水合版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AcH11")             # 纯碳酸插层型 AFm，水合度较高 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AcH9")              # 纯碳酸插层型 AFm，水合度较低 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AsH16")             # 硫酸型 AFm，结构为 C4ASH16 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH14")             # 硫酸插层 AFm 的低水合版本 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH12")             # C4ASH12，典型的单硫酸盐 AFm 相当于 monosulphate，常见于 OPC 与石膏反应后期
    # gemsk.supress_phase("C4AsH105")            # 水合程度更低的硫酸插层型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH9")              # 极低水合度硫酸型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    gemsk.supress_phase("Friedels")            # Friedel’s 盐，Cl⁻ 插层型 AFm，相当于 C4AClH10 Cl⁻ 存在时极易生成，尤其在含氯体系或使用海砂水泥中
    # gemsk.supress_phase("C2ASH55")             # 含硅 AFm 类型，C2ASH5.5，过渡相 低水合硅铝酸盐，硅含量较高或活性矿物体系中更常出现（如偏高岭土胶凝）
    # gemsk.supress_phase("C4FH13")              # 含铁型 AFm，Fe 替代部分 Al Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("C4Fc05H10")           # 混合铁铝型 AFm，水合度较低 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("C4FcH12")             # Fe 和 Al 共掺，水合度较高 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("monosulph-FeAl")      # 单硫酸盐型 AFm，含 Fe/Al 混合 Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    # gemsk.supress_phase("monosulph-AlFe")      # 同类结构，不同比例 Fe/Al  Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    # gemsk.supress_phase("SO4_OH_AFm")          # 硫酸根与 OH⁻ 共存插层型 AFm 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    # gemk.supress_phase("OH_SO4_AFm")          # OH⁻ 与硫酸根共存插层型 AFm，比例不同 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    # gemsk.supress_phase("C2AH75")              # 低水合度的 C2AH 相，常为过渡态 中间水合物，水化过程中的过渡相，尤其在实验初期较常见
    # gemsk.supress_phase("Kuzels")              # Kuzel’s 盐，Cl⁻ + SO₄²⁻ 双插层 AFm，特殊结构 Cl⁻+SO₄²⁻ 双插层型 AFm，较特殊，常见于海水或特殊氯盐体系

    
    # === Hydrogarnet（类水石榴石） ===
    # gemsk.supress_phase("C3AH6")               # 水合铝酸钙石榴石，C₃AH₆，稳定的水化产物之一 是 OPC 水化中铝酸钙类相的重要稳定产物（六方 → 立方结构转变）有在 高温（> 60°C） 条件下，水合铝酸盐类才会从 AFt（如 ettringite）/AFm 转变为 C₃AH₆。
    # gemlsk.supress_phase("C3FH6")               # 铁代铝的类石榴石结构矿物，C₃FH₆ Fe 掺杂型类石榴石，需要高 Fe₂O₃，普通 OPC 中不易形成
    # gemsk.supress_phase("C3FS0.84H4.32")       # 部分取代型含硅铁石榴石结构 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    # gemsk.supress_phase("C3(AF)S0.84H")        # 铁铝硅共掺的水化产物 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    # gemsk.supress_phase("C3FS1.34H3.32")       # 高硅含量的类石榴石矿物，结构更复杂 更复杂的高硅铁石榴石相，模拟中可能造成不必要相转移
    
    # === Al(OH)₃（铝氢氧化物） ===
    # gemsk.supress_phase("Al(OH)3am")           # 非晶态氢氧化铝，模拟初期水化或胶凝态存在 初期水化中可能有非晶铝羟胶出现，但难以定量，压不压看模拟目标
    # gemsk.supress_phase("Al(OH)3mic")          # 微晶氢氧化铝，具一定结晶度 微晶氢氧化铝，若系统中铝离子浓度较高，可析出；否则较少生成
    # gemsk.supress_phase("Gibbsite")            # 三水铝石，Al(OH)₃，自然界稳定相，常见于铝土矿 常见于自然铝土矿环境中，不是 OPC 水化产物，除非有偏高岭土或极强碱性体系
    
    # === Hydrotalcite（类水滑石） ===
    # gemsk.supress_phase("hydrotalc-pyro")      # 焙烧后再水化的水滑石（LDH）类矿物 焙烧型水滑石，需人工处理或Mg-Al来源，正常水泥中无来源
    # gemsk.supress_phase("OH-hydrotalcite")     # 氢氧根型水滑石，层状结构，具阴离子吸附功能 镁硅水合物，模拟低Ca胶凝体系时才出现（如 MgO + SiO₂ 体系），OPC 中基本不生成
    
    # === M-S-H（镁硅水合物） ===
    # gemsk.supress_phase("MSH")                 # 镁硅水合物，Mg-Si-H 凝胶态矿物，类似于低 Ca 的 C-S-H
    
    # === Zeolites（沸石类矿物） ===
    gemsk.supress_phase("ZeoliteP")            # 沸石 P，合成类沸石，常用于吸附/离子交换模拟 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gemsk.supress_phase("Natrolite")           # 钠沸石，Na₂Al₂Si₃O₁₀·2H₂O，自然沸石类型之一 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gemsk.supress_phase("Chabazite")           # 查巴茨沸石，具有较大孔径，阳离子交换能力强 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gemsk.supress_phase("ZeoliteX")            # 合成沸石 X，工业催化与分子筛常用 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gemsk.supress_phase("ZeoliteY")            # 沸石 Y，石化催化常用沸石，结构稳定 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。

    
    # === Clinkers（熟料矿物） ===
    # gemsk.supress_phase("Belite")         # 2CaO·SiO₂，低温稳定的硅酸二钙，相对水化慢
    # gemsk.supress_phase("Aluminate")      # 铝酸盐类，可能指 C₃A（三钙铝酸），对硫酸盐敏感
    # gemsk.supress_phase("Alite")          # 3CaO·SiO₂，高温形成，是水泥早期强度的主要来源
    # gemsk.supress_phase("Ferrite")        # 四钙铁铝酸盐，C₄AF，含铁熟料相
    gemsk.supress_phase("CA")             # 一钙铝酸盐，水化产物或熟料矿相 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gemsk.supress_phase("CA2")            # 二钙铝酸盐，水化速度较慢 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gemsk.supress_phase("lime")           # 氧化钙，游离石灰，可导致膨胀或不稳定  游离氧化钙（CaO），可能存在于部分高温水泥中，容易水化成 Portlandite。如果你不想模拟可能产生膨胀或碱反应，也可压制。
    gemsk.supress_phase("arcanite")       # K₂SO₄，硫酸钾，可能来自掺合料或污染源  这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gemsk.supress_phase("thenardite")     # Na₂SO₄，十水硫酸钠，高溶解度，可引起盐胀 这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gemsk.supress_phase("Na-oxide")       # 氧化钠，通常来自助熔剂或玻璃相成分 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    gemsk.supress_phase("K-oxide")        # 氧化钾，常出现在水泥熟料中，影响硅酸盐相稳定性 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    
    # === CaCO₃（碳酸钙类） ===
    gemsk.supress_phase("Aragonite")      # 霰石，CaCO₃ 的另一晶型，常在低温或有机环境中形成
    # gem.supress_phase("Calcite")        # 方解石，最稳定的 CaCO₃ 晶体，广泛存在于自然界 霰石为 CaCO₃ 的亚稳态晶型，需要低温或生物有机环境，OPC 水化中不易出现
    
    # === Portlandite（氢氧化钙） ===
    # gem.supress_phase("Portlandite")    # Ca(OH)₂，水泥水化的主要产物之一，提供碱性环境
    
    # === Sulfates（硫酸盐矿物） ===
    # gemk.supress_phase("Anhydrite")      # 无水硫酸钙，CaSO₄，常作为脱水后形式存在 若不使用脱水石膏作调凝剂可以压制，常见于特种水泥
    # gemk.supress_phase("Gypsum")         # 石膏，CaSO₄·2H₂O，常见调凝剂
    # gemk.supress_phase("hemihydrate")    # 半水石膏，CaSO₄·0.5H₂O，建筑石膏的主要形态 半水石膏主要是建筑石膏相，在 OPC 水泥中不是常见水化中间体

    
    # === Fe-phases（铁矿物） ===
    #原因：OPC 中 Fe₂O₃ 含量虽有，但远低于形成这些稳定铁矿相的需求。除非你使用钢渣、红泥、Fe 掺合料，否则这些相不会自然出现。
    gemsk.supress_phase("Iron")              # 金属铁，可能用于模拟原始还原性条件下的零价铁存在 
    gemsk.supress_phase("Fe-carbonate")      # 铁碳酸盐（可能为总称），例如菱铁矿等类矿物
    gemsk.supress_phase("Siderite")          # 菱铁矿，FeCO₃，典型的还原条件下形成的碳酸盐矿物
    gemsk.supress_phase("Hematite")          # 赤铁矿，Fe₂O₃，最稳定的三价铁氧化物，红褐色
    gemsk.supress_phase("Magnetite")         # 磁铁矿，Fe₃O₄，混合价态（Fe²⁺/Fe³⁺），具有强磁性
    gemsk.supress_phase("Ferrihydrite-am")   # 非晶铁水合氧化物，模拟新生成的纳米尺度氧化铁
    gemsk.supress_phase("Ferrihydrite-mc")   # 微晶型铁水合氧化物，Ferrihydrite 的结构化形式
    gemsk.supress_phase("Goethite")          # 针铁矿，α-FeOOH，常见于风化壤土和铁锈中
    gemsk.supress_phase("Pyrite")            # 黄铁矿，FeS₂，常见硫化物矿物，亮黄色金属光泽
    gemsk.supress_phase("Troilite")          # 斜方硫化亚铁，FeS，常见于陨石中，也存在于还原环境
    gemsk.supress_phase("Melanterite")       # 绿矾，FeSO₄·7H₂O，易溶于水，氧化还原环境中生成
    
    # === Mg-phases（镁矿物） ===
#OPC 中基本无 Mg 来源，除非有镁基掺合料（如 MgO、海砂、海水
    # gemk.supress_phase("Magnesite")         # 菱镁矿，MgCO₃，常见的镁碳酸盐矿物
    # gemk.supress_phase("Brucite")           # 氢氧化镁，Mg(OH)₂，在碱性环境中可析出，弱结晶性

    
    # === Mn-phases（锰矿物） ===
    gemsk.supress_phase("Rhodochrosite")       # 菱锰矿，MnCO₃，锰的主要碳酸盐矿物，粉红色晶体
    gemsk.supress_phase("Rhodochrosite-sy")    # 合成菱锰矿，模拟人工或非天然条件下形成的 MnCO₃
    gemsk.supress_phase("Hausmannite")         # 软锰矿，Mn₃O₄，混合价态氧化物，常见于风化环境中
    gemsk.supress_phase("Pyrolusite")          # 二氧化锰，MnO₂，天然最常见的锰氧化物，黑色矿物
    gemsk.supress_phase("Manganite")           # 水合氧化锰，MnO(OH)，一种含水氧化物，斜方晶系
    gemsk.supress_phase("Pyrochroite")         # 氢氧化锰，Mn(OH)₂，常见于还原环境或低氧条件下

    # === Other ===
    # === Clay minerals（粘土矿物） ===
    gemsk.supress_phase("Kaolinite")  # 高岭土，层状铝硅酸盐矿物，常见于土壤和沉积物中 常见于土壤或天然黏土系统，但在你的水泥+矿渣+粉煤灰+硅灰体系中不具形成条件。如果你没有添加天然黏土矿物，它不会生成
    
    # === Carbon phases（碳相） ===
    gemsk.supress_phase("Graphite")  # 石墨，碳的结晶形式，用于导电或润滑材料 石墨是纯碳的晶体，在任何正常水泥反应体系中都不会自然生成；除非模拟炭黑、还原环境，否则毫无意义
    
    # === Clinkers（熟料矿物） ===
    gemsk.supress_phase("Mayenite")  # 马叶矿，钙铝氧化物，在水泥熟料中形成于高温阶段 在高铝水泥或高温煅烧体系（如CSA或CA水泥）中可生成。你使用的是 OPC + 辅料体系，若未引入高铝熟料，则建议压制
    
    # === Carbonates（碳酸盐） ===
#白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gemsk.supress_phase("Dolomite-dis")  # 无序白云石，相较于有序白云石结构更松散 白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gemsk.supress_phase("Dolomite-ord")  # 有序白云石，CaMg(CO3)2，沉积岩中常见碳酸盐矿物 
    
    # === Sulfates（硫酸盐） ===
    gemsk.supress_phase("syngenite")  # 辛格石，K2Ca(SO4)2·H2O，一种钾钙复合硫酸盐矿物 K₂Ca(SO₄)₂·H₂O，特殊的钾钙硫酸盐，需高浓度 K⁺ + SO₄²⁻ 才可能生成，不属于正常 OPC 水化产物，可能带来模型不稳定或假阳离子迁移结果
    
    # === Elemental phases（元素相） ===
    # gemk.supress_phase("Sulphur")  # 硫，常以单质形式出现于火山或矿物环境中 单质硫只在极强还原或火山环境中形成，水泥体系无形成机制
    
    # === Silica phases（硅相） ===
    # gemk.supress_phase("Quartz")  # 石英，最常见的二氧化硅晶体结构 粉煤灰和矿渣中可能含有惰性石英，但通常不参与反应。若只关注活性反应产物，可压制；若要模拟总固相矿物平衡，可保留
    # gem.supress_phase("Silica-amorph")  # 非晶态二氧化硅，模拟水化产物或凝胶硅质成分 是硅灰、粉煤灰等活性成分形成的主要凝胶态结构之一，对 C-S-H（或 C-A-S-H）模拟至关重要，强烈建议保留
    
    # === Titanium oxides（钛氧化物） ===
#Ti 元素在你体系中基本不存在，TiO₂ 常用于颜料或光催化水泥，与结构矿物无关，压制无影响
    gemsk.supress_phase("Ti(alpha)")  # α-钛，六方结构的金属钛，在氧化之前的中间物相
    gemsk.supress_phase("TiO2(am_hyd)")  # 非晶水合二氧化钛，常见于纳米材料和胶体系统中

    
    gemsk.add_species_amt("O2@",0.001, units = 'kg')#氧化所以总质量为多余0.1g
    equilibration = gemsk.equilibrate()
    if equilibration[0] == 'F' or equilibration[0] == 'B':
        if fail:
            print('--- Equilibration Failed ---')
            print('--- Trying Again with more H2O ---')

            # 添加额外水量后再试一次
            gemsk.add_species_amt("H2O@", 0.001, units="kg")
            equilibration = gemsk.equilibrate()

            # 再次失败，打印问题配方
            if equilibration[0] == 'F' or equilibration[0] == 'B':
                print(' **Fail to equilibrate for recipe:**')
                print(recipe)  # 打印失败的配方
                return None  # 继续运行，而不是中断
            else:
                print('--- Equilibration OK after adding H2O ---')
                return process_gems_result(gemsk)

        else:
            print('**Equilibration Failed, but continuing...**')
            print(recipe)  # 直接打印失败的配方
            return {},{},{},{},{},{}, {},{},{},{},{} # 默认值为空字典
  # 不中断程序

    else:
        return process_gems_result(gemsk)
    
def alpha(t: float, SCM_type):
    ''' Compute the alpha in 5PL model in https://cemgems.org/tutorial/enhanced/redefining-processes/#b-the-5pl-4pl-logistic-function-model
        Used for the hydration of SCM (only implemented for GGBFS, Calcined Clay and Fly Ash)
    '''
    if t == 0.:
        alpha= dic_5lp[SCM_type]['A']
    else:
        constant = dic_5lp[SCM_type]
        alpha = constant['D'] + (constant['A']-constant['D'])/(1+(t/constant['C'])**constant['B'])**constant['G']
    return alpha




def hydration(recipe, output_times, fail=False):
    pk = parrot_killoh(recipe['wc'], recipe['RH'], recipe['T'], recipe['fineness'])
    
    input_file = 'gems_files/MySystem2-dat.lst'
    gemsk = GEMS(input_file)
     # === AFt-phases（三钙铝石类，AFt） ===
    # gemsk.supress_phase("ettringite-AlFe")     # 含铝铁的钙矾石，AFt 结构变体之一 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    # gemsk.supress_phase("ettringite-FeAl")     # 含铁铝的钙矾石，另一种AFt结构混合型 需要较高 Fe³⁺ 参与，OPC 中 Fe 含量较低，仅在富铁体系或引入钢渣/红泥时可能形成
    # gemsk.supress_phase("ettringite")          # 钙矾石，C6AŜ3H32，典型的硫铝酸盐水化产物
    # # gemsk.supress_phase("SO4_CO3_AFt")         # 硫酸根-碳酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    # gemsk.supress_phase("CO3_SO4_AFt")         # 碳酸根-硫酸根共存的 AFt 结构变体 需碳酸根参与（如CO₂充足、碳酸盐污染），正常 OPC 水化中不常见，属于混合阴离子 AFt
    # gemsk.supress_phase("C6AsH13")             # 硫酸型AFt矿物变体，可能水合度不同 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    # gemsk.supress_phase("C6AsH9")              # 更低水合度的硫酸型AFt矿物 属于理论/推导型 AFt 相（低水合度钙矾石变体），实验中罕有报道，不常自然生成
    gemsk.supress_phase("thaumasite")          # 硅钙矾石（Thaumasite），Si 替代 SO₄ 形成的 AFt 型结构 需低温（<15°C）+高 CO₂ + Si源（如硅灰）+ SO₄²⁻，在 OPC 正常条件下极少生成，但在寒冷潮湿环境下可能导致腐蚀问题（称为“thaumasite 硫酸盐攻击”）
    
    # === AFm-Silica（硅铝酸盐型 AFm） ===
    # gemsk.supress_phase("straetlingite")       # 含硅 AFm，相当于 C2ASH8，常见于硅含量较高体系 在硅铝材料体系（如含火山灰、偏高岭土的胶凝材料）中可见，在纯 OPC 中较少见，但并非罕见
    
    # === AFm-phases（单盐型水化铝酸盐 AFm） ===
    # gemsk.supress_phase("C4AH19")              # 高水合度铝酸钙水合物，早期水化产物 高水合度，早期生成，
    # gemsk.supress_phase("C4AH13")              # 中等水合度 AFm，常在水化过程中生成 ，是常见稳定产物之一
    # gemsk.supress_phase("C4AH11")              # 较低水合度的 AFm，相对更稳定 ，相对稳定，晚期可能形成
    # gemsk.supress_phase("CAH10")               # 一种低温条件下形成的 AFm 类型 低温环境或高水灰比条件下更易生成，实验中难观察
    # gemsk.supress_phase("C4Ac0.5H12")          # 碳酸-氢氧共插层型 AFm，相对较高水合度（CO₃²⁻ + OH⁻）需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4Ac0.5H105")         # 同上，略低水合度版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4Ac0.5H9")           # 碳酸型 AFm 的低水合版本 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AcH11")             # 纯碳酸插层型 AFm，水合度较高 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AcH9")              # 纯碳酸插层型 AFm，水合度较低 需要大量 CO₃²⁻（碳酸根），如暴露在 CO₂ 浓环境中或碳酸盐掺合料时才会形成
    # gemsk.supress_phase("C4AsH16")             # 硫酸型 AFm，结构为 C4ASH16 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH14")             # 硫酸插层 AFm 的低水合版本 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH12")             # C4ASH12，典型的单硫酸盐 AFm 相当于 monosulphate，常见于 OPC 与石膏反应后期
    # gemsk.supress_phase("C4AsH105")            # 水合程度更低的硫酸插层型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    # gemsk.supress_phase("C4AsH9")              # 极低水合度硫酸型 AFm 理论存在的不同水合度 AFm，实验中不常精确分辨，仅模拟细化时保留
    gemsk.supress_phase("Friedels")            # Friedel’s 盐，Cl⁻ 插层型 AFm，相当于 C4AClH10 Cl⁻ 存在时极易生成，尤其在含氯体系或使用海砂水泥中
    # gemsk.supress_phase("C2ASH55")             # 含硅 AFm 类型，C2ASH5.5，过渡相 低水合硅铝酸盐，硅含量较高或活性矿物体系中更常出现（如偏高岭土胶凝）
    # gemsk.supress_phase("C4FH13")              # 含铁型 AFm，Fe 替代部分 Al Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("C4Fc05H10")           # 混合铁铝型 AFm，水合度较低 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("C4FcH12")             # Fe 和 Al 共掺，水合度较高 Fe 掺杂型 AFm，需要显著的 Fe₂O₃ 来源，如钢渣等
    # gemsk.supress_phase("monosulph-FeAl")      # 单硫酸盐型 AFm，含 Fe/Al 混合 Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    # gemsk.supress_phase("monosulph-AlFe")      # 同类结构，不同比例 Fe/Al  Fe/Al 掺杂 AFm，需高 Fe 含量，通常 OPC 无法提供
    # gemsk.supress_phase("SO4_OH_AFm")          # 硫酸根与 OH⁻ 共存插层型 AFm 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    # gemk.supress_phase("OH_SO4_AFm")          # OH⁻ 与硫酸根共存插层型 AFm，比例不同 混合阴离子型 AFm，现实中易出现于硫酸根含量不足或后期平衡过程中
    # gemsk.supress_phase("C2AH75")              # 低水合度的 C2AH 相，常为过渡态 中间水合物，水化过程中的过渡相，尤其在实验初期较常见
    # gemsk.supress_phase("Kuzels")              # Kuzel’s 盐，Cl⁻ + SO₄²⁻ 双插层 AFm，特殊结构 Cl⁻+SO₄²⁻ 双插层型 AFm，较特殊，常见于海水或特殊氯盐体系

    
    # === Hydrogarnet（类水石榴石） ===
    # gemsk.supress_phase("C3AH6")               # 水合铝酸钙石榴石，C₃AH₆，稳定的水化产物之一 是 OPC 水化中铝酸钙类相的重要稳定产物（六方 → 立方结构转变）有在 高温（> 60°C） 条件下，水合铝酸盐类才会从 AFt（如 ettringite）/AFm 转变为 C₃AH₆。
    # gemlsk.supress_phase("C3FH6")               # 铁代铝的类石榴石结构矿物，C₃FH₆ Fe 掺杂型类石榴石，需要高 Fe₂O₃，普通 OPC 中不易形成
    # gemsk.supress_phase("C3FS0.84H4.32")       # 部分取代型含硅铁石榴石结构 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    # gemsk.supress_phase("C3(AF)S0.84H")        # 铁铝硅共掺的水化产物 含 Fe 和 Si 的类石榴石，需特殊 Fe/Si 环境
    # gemsk.supress_phase("C3FS1.34H3.32")       # 高硅含量的类石榴石矿物，结构更复杂 更复杂的高硅铁石榴石相，模拟中可能造成不必要相转移
    
    # === Al(OH)₃（铝氢氧化物） ===
    # gemsk.supress_phase("Al(OH)3am")           # 非晶态氢氧化铝，模拟初期水化或胶凝态存在 初期水化中可能有非晶铝羟胶出现，但难以定量，压不压看模拟目标
    # gemsk.supress_phase("Al(OH)3mic")          # 微晶氢氧化铝，具一定结晶度 微晶氢氧化铝，若系统中铝离子浓度较高，可析出；否则较少生成
    # gemsk.supress_phase("Gibbsite")            # 三水铝石，Al(OH)₃，自然界稳定相，常见于铝土矿 常见于自然铝土矿环境中，不是 OPC 水化产物，除非有偏高岭土或极强碱性体系
    
    # === Hydrotalcite（类水滑石） ===
    # gemsk.supress_phase("hydrotalc-pyro")      # 焙烧后再水化的水滑石（LDH）类矿物 焙烧型水滑石，需人工处理或Mg-Al来源，正常水泥中无来源
    # gemsk.supress_phase("OH-hydrotalcite")     # 氢氧根型水滑石，层状结构，具阴离子吸附功能 镁硅水合物，模拟低Ca胶凝体系时才出现（如 MgO + SiO₂ 体系），OPC 中基本不生成
    
    # === M-S-H（镁硅水合物） ===
    # gemsk.supress_phase("MSH")                 # 镁硅水合物，Mg-Si-H 凝胶态矿物，类似于低 Ca 的 C-S-H
    
    # === Zeolites（沸石类矿物） ===
    gemsk.supress_phase("ZeoliteP")            # 沸石 P，合成类沸石，常用于吸附/离子交换模拟 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gemsk.supress_phase("Natrolite")           # 钠沸石，Na₂Al₂Si₃O₁₀·2H₂O，自然沸石类型之一 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gemsk.supress_phase("Chabazite")           # 查巴茨沸石，具有较大孔径，阳离子交换能力强 天然沸石，需长时间与硅铝胶体反应，或特定温压条件（如高温水热养护），常出现在水泥/矿渣长期碳化后或沸石水泥中，OPC 水化中基本没
    gemsk.supress_phase("ZeoliteX")            # 合成沸石 X，工业催化与分子筛常用 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。
    gemsk.supress_phase("ZeoliteY")            # 沸石 Y，石化催化常用沸石，结构稳定 人工合成型沸石，用于催化和吸附模拟，不是水泥中自然形成相。

    
    # === Clinkers（熟料矿物） ===
    # gemsk.supress_phase("Belite")         # 2CaO·SiO₂，低温稳定的硅酸二钙，相对水化慢
    # gemsk.supress_phase("Aluminate")      # 铝酸盐类，可能指 C₃A（三钙铝酸），对硫酸盐敏感
    # gemsk.supress_phase("Alite")          # 3CaO·SiO₂，高温形成，是水泥早期强度的主要来源
    # gemsk.supress_phase("Ferrite")        # 四钙铁铝酸盐，C₄AF，含铁熟料相
    gemsk.supress_phase("CA")             # 一钙铝酸盐，水化产物或熟料矿相 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gemsk.supress_phase("CA2")            # 二钙铝酸盐，水化速度较慢 属于高铝水泥（CAC）体系的熟料相，不是普通 OPC 的常规成分。
    gemsk.supress_phase("lime")           # 氧化钙，游离石灰，可导致膨胀或不稳定  游离氧化钙（CaO），可能存在于部分高温水泥中，容易水化成 Portlandite。如果你不想模拟可能产生膨胀或碱反应，也可压制。
    gemsk.supress_phase("arcanite")       # K₂SO₄，硫酸钾，可能来自掺合料或污染源  这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gemsk.supress_phase("thenardite")     # Na₂SO₄，十水硫酸钠，高溶解度，可引起盐胀 这两个是可溶盐类（K₂SO₄ 和 Na₂SO₄），容易造成“盐胀”或不合理离子迁移，通常在水化模型中被人为排除。
    gemsk.supress_phase("Na-oxide")       # 氧化钠，通常来自助熔剂或玻璃相成分 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    gemsk.supress_phase("K-oxide")        # 氧化钾，常出现在水泥熟料中，影响硅酸盐相稳定性 模拟助熔剂来源时会用，但不是矿物相，通常也不会出现在平衡产物中，建议压制。
    
    # === CaCO₃（碳酸钙类） ===
    gemsk.supress_phase("Aragonite")      # 霰石，CaCO₃ 的另一晶型，常在低温或有机环境中形成
    # gem.supress_phase("Calcite")        # 方解石，最稳定的 CaCO₃ 晶体，广泛存在于自然界 霰石为 CaCO₃ 的亚稳态晶型，需要低温或生物有机环境，OPC 水化中不易出现
    
    # === Portlandite（氢氧化钙） ===
    # gem.supress_phase("Portlandite")    # Ca(OH)₂，水泥水化的主要产物之一，提供碱性环境
    
    # === Sulfates（硫酸盐矿物） ===
    # gemk.supress_phase("Anhydrite")      # 无水硫酸钙，CaSO₄，常作为脱水后形式存在 若不使用脱水石膏作调凝剂可以压制，常见于特种水泥
    # gemk.supress_phase("Gypsum")         # 石膏，CaSO₄·2H₂O，常见调凝剂
    # gemk.supress_phase("hemihydrate")    # 半水石膏，CaSO₄·0.5H₂O，建筑石膏的主要形态 半水石膏主要是建筑石膏相，在 OPC 水泥中不是常见水化中间体

    
    # === Fe-phases（铁矿物） ===
    #原因：OPC 中 Fe₂O₃ 含量虽有，但远低于形成这些稳定铁矿相的需求。除非你使用钢渣、红泥、Fe 掺合料，否则这些相不会自然出现。
    gemsk.supress_phase("Iron")              # 金属铁，可能用于模拟原始还原性条件下的零价铁存在 
    gemsk.supress_phase("Fe-carbonate")      # 铁碳酸盐（可能为总称），例如菱铁矿等类矿物
    gemsk.supress_phase("Siderite")          # 菱铁矿，FeCO₃，典型的还原条件下形成的碳酸盐矿物
    gemsk.supress_phase("Hematite")          # 赤铁矿，Fe₂O₃，最稳定的三价铁氧化物，红褐色
    gemsk.supress_phase("Magnetite")         # 磁铁矿，Fe₃O₄，混合价态（Fe²⁺/Fe³⁺），具有强磁性
    gemsk.supress_phase("Ferrihydrite-am")   # 非晶铁水合氧化物，模拟新生成的纳米尺度氧化铁
    gemsk.supress_phase("Ferrihydrite-mc")   # 微晶型铁水合氧化物，Ferrihydrite 的结构化形式
    gemsk.supress_phase("Goethite")          # 针铁矿，α-FeOOH，常见于风化壤土和铁锈中
    gemsk.supress_phase("Pyrite")            # 黄铁矿，FeS₂，常见硫化物矿物，亮黄色金属光泽
    gemsk.supress_phase("Troilite")          # 斜方硫化亚铁，FeS，常见于陨石中，也存在于还原环境
    gemsk.supress_phase("Melanterite")       # 绿矾，FeSO₄·7H₂O，易溶于水，氧化还原环境中生成
    
    # === Mg-phases（镁矿物） ===
#OPC 中基本无 Mg 来源，除非有镁基掺合料（如 MgO、海砂、海水
    # gemk.supress_phase("Magnesite")         # 菱镁矿，MgCO₃，常见的镁碳酸盐矿物
    # gemk.supress_phase("Brucite")           # 氢氧化镁，Mg(OH)₂，在碱性环境中可析出，弱结晶性

    
    # === Mn-phases（锰矿物） ===
    gemsk.supress_phase("Rhodochrosite")       # 菱锰矿，MnCO₃，锰的主要碳酸盐矿物，粉红色晶体
    gemsk.supress_phase("Rhodochrosite-sy")    # 合成菱锰矿，模拟人工或非天然条件下形成的 MnCO₃
    gemsk.supress_phase("Hausmannite")         # 软锰矿，Mn₃O₄，混合价态氧化物，常见于风化环境中
    gemsk.supress_phase("Pyrolusite")          # 二氧化锰，MnO₂，天然最常见的锰氧化物，黑色矿物
    gemsk.supress_phase("Manganite")           # 水合氧化锰，MnO(OH)，一种含水氧化物，斜方晶系
    gemsk.supress_phase("Pyrochroite")         # 氢氧化锰，Mn(OH)₂，常见于还原环境或低氧条件下

    # === Other ===
    # === Clay minerals（粘土矿物） ===
    gemsk.supress_phase("Kaolinite")  # 高岭土，层状铝硅酸盐矿物，常见于土壤和沉积物中 常见于土壤或天然黏土系统，但在你的水泥+矿渣+粉煤灰+硅灰体系中不具形成条件。如果你没有添加天然黏土矿物，它不会生成
    
    # === Carbon phases（碳相） ===
    gemsk.supress_phase("Graphite")  # 石墨，碳的结晶形式，用于导电或润滑材料 石墨是纯碳的晶体，在任何正常水泥反应体系中都不会自然生成；除非模拟炭黑、还原环境，否则毫无意义
    
    # === Clinkers（熟料矿物） ===
    gemsk.supress_phase("Mayenite")  # 马叶矿，钙铝氧化物，在水泥熟料中形成于高温阶段 在高铝水泥或高温煅烧体系（如CSA或CA水泥）中可生成。你使用的是 OPC + 辅料体系，若未引入高铝熟料，则建议压制
    
    # === Carbonates（碳酸盐） ===
#白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gemsk.supress_phase("Dolomite-dis")  # 无序白云石，相较于有序白云石结构更松散 白云石是钙镁双盐类碳酸盐，需大量 Mg²⁺ 和 CO₃²⁻，在你体系中无 Mg 来源，极不易形成
    gemsk.supress_phase("Dolomite-ord")  # 有序白云石，CaMg(CO3)2，沉积岩中常见碳酸盐矿物 
    
    # === Sulfates（硫酸盐） ===
    gemsk.supress_phase("syngenite")  # 辛格石，K2Ca(SO4)2·H2O，一种钾钙复合硫酸盐矿物 K₂Ca(SO₄)₂·H₂O，特殊的钾钙硫酸盐，需高浓度 K⁺ + SO₄²⁻ 才可能生成，不属于正常 OPC 水化产物，可能带来模型不稳定或假阳离子迁移结果
    
    # === Elemental phases（元素相） ===
    # gemk.supress_phase("Sulphur")  # 硫，常以单质形式出现于火山或矿物环境中 单质硫只在极强还原或火山环境中形成，水泥体系无形成机制
    
    # === Silica phases（硅相） ===
    # gemsk.supress_phase("Quartz")  # 石英，最常见的二氧化硅晶体结构 粉煤灰和矿渣中可能含有惰性石英，但通常不参与反应。若只关注活性反应产物，可压制；若要模拟总固相矿物平衡，可保留
    # gemsk.supress_phase("Silica-amorph")  # 非晶态二氧化硅，模拟水化产物或凝胶硅质成分 是硅灰、粉煤灰等活性成分形成的主要凝胶态结构之一，对 C-S-H（或 C-A-S-H）模拟至关重要，强烈建议保留
    
    # === Titanium oxides（钛氧化物） ===
#Ti 元素在你体系中基本不存在，TiO₂ 常用于颜料或光催化水泥，与结构矿物无关，压制无影响
    gemsk.supress_phase("Ti(alpha)")  # α-钛，六方结构的金属钛，在氧化之前的中间物相
    gemsk.supress_phase("TiO2(am_hyd)")  # 非晶水合二氧化钛，常见于纳米材料和胶体系统中
    
    gems_vol = {}
    gems_vol_frac = {}
    gems_phase_masses = {}
    remaining_vol_dict = {}
    remaining_mass_dict = {}
    gems_pH_dic = {}
    density = []
    element_CSHQ = {}
    element_aq = {}
    aq_composition = {}
    gems_B = {}
    gems_CSH_species_mass = {}
    gems_CSH_species_volumes = {}
    intial_vol = {}
    index_to_drop = []

    for i in range(len(output_times)):
        fineness = recipe['fineness']
        CO2 = recipe['CO2']
        wc = recipe['wc']
        RH = recipe['RH']
        T = recipe['T']
        pk = parrot_killoh(wc, RH, T, fineness)
        
        clink_phases = {key: recipe[key] for key in ['C3S', 'C2S', 'C3A', 'C4AF']}
        CSH2 = recipe['CSH2']
        all_species = clink_phases.copy()
        all_species["H2O@"] = wc * 100
        all_species["CO2"] = CO2
        all_species["Gp"] = CSH2
        
        for name in all_species.keys():
            all_species[name] *= 1e-3
        
        gemsk.T = T + 273.15
        gemsk.add_multiple_species_amt(all_species, units="kg")
        gemsk.add_species_amt("O2@", 0.001, units='kg')

        if i > 0:
            gemsk.warm_start()

        for phase in clink_phases:
            gemsk.species_lower_bound(phase, clink_phases[phase] * (1 - pk[phase][i]) * 1e-3, units="kg")

        for component in recipe.keys():
            # if component == 'limestone':
            #     continue
            if component not in ['C3S', 'C2S', 'C3A', 'C4AF', 'wc', 'CSH2', 'T', 'RH', 'fineness', 'CO2']:
                try:
                    if component in SCM.keys():
                        value = recipe[component] * 1e-3 * alpha(output_times[i], component) * 1e-2
                    else:
                        value = recipe[component] * 1e-3
                    
                    remaining_mass = (recipe[component] * 1e-3) - value
        
                    if component in SCM_density.keys():
                        remaining_vol = remaining_mass / SCM_density[component]
                    else:
                        remaining_vol = remaining_mass / 2000
        
                    if output_times[i] not in remaining_vol_dict:
                        remaining_vol_dict[output_times[i]] = {}
                        remaining_mass_dict[output_times[i]] = {}
                        
                    remaining_vol_dict[output_times[i]][component] = remaining_vol
                    remaining_mass_dict[output_times[i]][component] = remaining_mass
        
                except KeyError:
                    value = recipe[component] * 1e-3

                if SCM[component][1] == 'mol':
                    gemsk.add_amt_from_formula(SCM[component][0], value, units='mol')
                elif SCM[component][1] == 'kg':
                    gemsk.add_amt_from_formula(SCM[component][0], value, units='kg')
                else:
                    raise UnitsError('the unit is not taken into account please select "mol" or "kg"')

        equilibration = gemsk.equilibrate()
        
        if equilibration[0] in ['F', 'B']:
            index_to_drop.append(i)
            # print(f"Equilibration failed for time step {output_times[i]}")
            return index_to_drop, None, None, None, None, None, None, None, None, None
        else:
            gems_phase_masses[output_times[i]] = gemsk.phase_masses
            gems_CSH_species_mass[output_times[i]] = gemsk.cshq_species_masses
            gems_CSH_species_volumes[output_times[i]] = gemsk.cshq_species_volumes
            gems_vol[output_times[i]] = gemsk.phase_volumes
            gems_pH_dic[output_times[i]] = gemsk.pH
            gems_B[output_times[i]] = gemsk.b.copy()
            gemsk.clear()

            oh_value = gemsk.aq_composition.get('OH-', np.float64(0.0))
            element_aq_list = gemsk.aq_element.tolist()
            element_aq_list.append(oh_value)
            element_aq[output_times[i]] = np.array(element_aq_list)
            element_CSHQ[output_times[i]] = gemsk.reset_CSHQ_composition
            intial_vol[output_times[i]] = gemsk.system_volume

            if output_times[i] in remaining_vol_dict:
                gems_vol[output_times[i]] = {**remaining_vol_dict[output_times[i]], **gems_vol[output_times[i]]}
            if output_times[i] in remaining_mass_dict:
                gems_phase_masses[output_times[i]] = {**remaining_mass_dict[output_times[i]], **gems_phase_masses[output_times[i]]}

            for key in gems_vol[output_times[i]].keys():
                gems_vol[output_times[i]][key] /= 1
            for key in gems_phase_masses[output_times[i]].keys():
                gems_phase_masses[output_times[i]][key] /= 1
            for key in gems_CSH_species_mass[output_times[i]].keys():
                gems_CSH_species_mass[output_times[i]][key] /= 1
            for key in gems_CSH_species_volumes[output_times[i]].keys():
                gems_CSH_species_volumes[output_times[i]][key] /= 1
            for key in gems_phase_masses[output_times[i]].keys():
                gems_phase_masses[output_times[i]][key] /= 0.001
                density.append(gemsk.system_mass / gemsk.system_volume)
            # print( gems_vol)

    return index_to_drop, gems_vol, gems_phase_masses, density, gems_pH_dic, gems_B, element_CSHQ, gems_CSH_species_mass, gems_CSH_species_volumes, element_aq


def composition_only(recipe, output_times):
    '''计算体系内的元素量（gems.b），并考虑 SCM 在不同时间的反应程度（alpha）'''
    
    input_file = 'gems_files/MySystem2-dat.lst'
    gemsk = GEMS(input_file)
    
    all_B = {}
    
    for i, time in enumerate(output_times):
        gemsk.clear()
        selected_keys = ['C3S', 'C2S', 'C3A', 'C4AF', 
                     'limestone', 'silica_fume', 'GGBFS', 'fly_ash', 'calcined_clay',
                     'CSH2',]
        # H2O 和 CSH2 特殊处理
        water_kg = recipe['wc'] * 0 * 1e-3
        CSH2_kg = recipe['CSH2'] * 1e-3
        gemsk.add_species_amt("H2O@", water_kg, units="kg")
        gemsk.add_species_amt("Gp", CSH2_kg, units="kg")
        
        # 添加熟料相和 SCM
        for key in selected_keys:
            if key in ['CSH2']:
                continue  # 已添加
            if key in recipe:
                mass_kg = recipe[key] * 1e-3
                
                if key in dic_5lp:
                    alpha_val = alpha(time, key) * 1e-2
                    mass_kg *= alpha_val  # 考虑反应程度

                if key in ['C3S', 'C2S', 'C3A', 'C4AF']:
                    gemsk.add_species_amt(key, mass_kg, units='kg')
                else:
                    gemsk.add_amt_from_formula(SCM[key][0], mass_kg, units='kg')
        
        # gemsk.add_species_amt("O2@", 1e-6, units="kg")  # 减少红氧问题
        
        B = gemsk.b.copy()
        all_B[time] = dict(zip(element_names[:-1], B))  # 去掉 OH-，因为没被加
        
    return all_B


