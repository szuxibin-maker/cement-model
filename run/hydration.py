import numpy as np
from scipy.integrate import solve_ivp
from run.GEMSCalc import GEMS
import matplotlib.pylab as plt
import cmocean
#from bqplot import pyplot as plt
import matplotlib.patches as mpatches
from run import Time_points
# # dictionary for the recipe of some SCM from CEMGEMS CEM-II-BV based on 100 g of each SCM  always give formula in moles!!!
# SCM = { 'limestone': ({'Al' : 0.052961280988, 'C' : 0.85663128045, 'Ca': 0.75431457236,'Fe' : 0.025048813876, 'K': 0.012739394454,'Mg' : 0.044660136362, 'Na': 0.016134497168, 'O': 3.0564427725, 'Si':0.20637670739},  'mol'),
#         'silica_fume':({'Al' : 0.005883425509, 'Ca': 0.002674346509,'Fe' : 0.0002504387073, 'K': 0.0063684402277,'Mg' :0.00099225162461, 'Mn': 0.0001266581442,'Na' : 0.0029036363721, 'O': 3.3128876396,'P' : 0.00028174463078, 'Si': 1.6423635207,'Ti' : 0.0048814082984}, 'mol'),
#        'GGBFS' : ({'Al' : 0.25068590354, 'C' : 0.026137740473, 'Ca': 0.7400381436, 'Fe': 0.0062268216656, 'H': 0.24594931045, 'K': 0.0066234897348, 'Mg':0.18043867977, 'Na': 0.0050333680703, 'O':2.568871619, 'S': 0.034432524341, 'Si': 0.54806878508, 'Ti':0.0064682389664},'mol'), 
#        'calcined_clay':({'Al' : 0.86001188208, 'Ca' : 0.0017850310446, 'Fe': 0.0037610756424, 'H': 0.16669147874, 'K': 0.0021253535161, 'Na':0.009690369309, 'O': 3.1677258055, 'P':0.0028208175955, 'S':0.0012502298115, 'Si':0.86631529281, 'Ti': 0.018797209003},'mol'), 
#        'fly_ash': ({'Al' : 0.35503673551, 'C' : 0.016132843743, 'Ca': 0.17672003338, 'Fe': 0.11522454383, 'K': 0.051594547539, 'Mg':0.087831601512, 'Na': 0.03259168428, 'O':2.8800652865,'S': 0.0049959283184, 'Si':  0.91038757213},'mol'),
#       'Pozzolan': ({'Al' : 0.3500773, 'Ca' : 0.1609977,  'Fe': 0.05360106, 'H':0.3377904, 'K': 0.1611911, 'Mg':0.03044469, 'Na': 0.09915054, 'O':2.902273,'S': 0.008098946, 'Si':  0.8813103},'mol'),
#        'biochar': ({'Al' : 0.01473923, 'Ca' : 1.076925,  'Fe': 0.03687062,  'K': 0.4672641, 'Mg':0.1200335, 'Na': 0.08669282, 'O':1.781338,'P':0.02683121,'S': 0.04139704, 'Si':  0.01935863},'mol')
#       }
# Carbonation = {'Co2':({'C':2.2727273, 'O':4.5454545 },'mol')}

# '''
# Formula in %mass for 100g: 
# - Limestone: CaO: 42.3, Al2O3: 2.7, Fe2O3: 2, SiO2: 12.4, Na2O: 0.5, MgO: 1.8, K2O: 0.6, CO2: 37.7

# - Silica Fume: CaO: 0.14997, MgO: 0.039992, K2O: 0.29994, Na2O: 0.089982, TiO2: 0.38992, Mn2O3: 0.009998, P2O5: 0.019996, Al2O3: 0.29994, Fe2O3: 0.019996, SiO2: 98.68

# - GGBFS: SiO2: 32.93, Al2O3: 12.78, Fe2O3: 0.49717, CaO: 41.499, MgO: 7.2724, SO3: 0.14623, Na2O: 0.15598, K2O: 0.31195, H2S: 1.1113, H2O: 1.628, TiO2: 0.51667, CO2: 1.1503

# - calcined clay: SiO2: 52.052, Al2O3: 43.844, Fe2O3: 0.3003, CaO: 0.1001, SO3: 0.1001, Na2O: 0.3003, K2O: 0.1001, TiO2: 1.5015, P2O5: 0.2002, H2O: 1.5015 

# - fly ash: Ca0: 9.91, Al2O3: 18.1, Fe2O2: 9.2, SiO2: 54.7, CO2: 0.71, MgO: 3.54, K2O: 2.43, Na2O: 1.01, SO2: 0.4, P2O5: 0, C: 0
# -Pozzolan ：CaO: 9.05, Al2O3: 17.89, Fe2O3: 4.29, SiO2: 53.08, H20:3.05 , MgO: 1.23, K2O: 7.61, Na2O: 3.08, SO3: 0.65, TiO2: 0.31
# - biochar (forest wood): CaO: 58.67, Al2O3: 0.73, Fe2O3: 2.86, SiO2: 1.13, MgO: 4.7, K2O: 21.38, Na2O: 2.61, SO3: 3.22,Cl:0.47,P2O5: 1.85
#  reference: Effect of cement partial substitution by waste-based biochar in mortars properties
# '''



# # dictionary for the 5PL model (values found in cemgem)
# dic_5lp = {'GGBFS' : {'A':0, 'B':1, 'C':7, 'D':70, 'G':1},
#            'calcined_clay' : {'A':0, 'B':1, 'C':5, 'D':50, 'G':1},
#            'fly_ash' : {'A':0, 'B':0.7, 'C':85.1, 'D':60, 'G':1},
#           }

# density = {'GGBFS' : {2.88},
#            'calcined_clay' : {2.3},
#            'fly_ash' : {2.49},
#           'Silica_fume' : {2.2},
#           'Pozzolana ' : {2},
#           'Limestone':{2.709888},}



def bouges_composition(CaO,SiO2,Al2O3,Fe2O3,SO3):
    """
    computation of clinker phases from oxides using bouge's conversion
    """
    clink_phases={}
    clink_phases["C3S"]=C3S=4.07*CaO -7.6*SiO2-6.72*Al2O3-1.43*Fe2O3-2.85*SO3 
    clink_phases["C2S"]=C2S=2.87*SiO2-0.75*C3S 
    clink_phases["C3A"]=C3A=2.65*Al2O3-1.69*Fe2O3 
    clink_phases["C4AF"]=C4AF=3.04*Fe2O3
    return clink_phases

#parameters for parrot and killoh model source: lothenbach et al. 2008
K1 = {"C3S":1.5,"C2S":0.5,"C3A":1.,"C4AF":0.37}
N1 = {"C3S":0.7,"C2S":1.0,"C3A":0.85,"C4AF":0.7}
K2 = {"C3S":0.05,"C2S":0.02,"C3A":0.04,"C4AF":0.015}
K3 = {"C3S":1.1,"C2S":0.7,"C3A":1.0,"C4AF":0.4}
N3 = {"C3S":3.3,"C2S":5.0,"C3A":3.2,"C4AF":3.7}
H  = {"C3S":1.8,"C2S":1.35,"C3A":1.6,"C4AF":1.45}
Ea = {"C3S":41570,"C2S":20785,"C3A":54040,"C4AF":34087} #from barbara's excel sheets in paper only 2 leading digits are given
T0 = 25#c
ref_fineness = 385

# output_times=[0., 7, 14, 28, 90] #range for output of time [days]. Always at 0 for initial conditions
output_times=Time_points.output_times

"""
This cell provides implementation of parrot and killoh rate model implemented above
"""
def nuclandgrowth(alpha_t,K1,N1,fineness,ref_fineness):
    """
    Nucleation and growth rate equation of parrot and killoh model
    """
    ln = np.log
    rt = (K1/N1) * (1-alpha_t) * (-ln(1-alpha_t))**(1-N1) * (fineness/ref_fineness)
    return rt

def diffusion(alpha_t,K2):
    """
    Diffusion rate equation of parrot and killoh model
    """
    eps = 1e-12
    alpha_safe = np.clip(alpha_t, eps, 1 - eps)
    rt = K2 * (1 - alpha_safe)**(2/3) / (alpha_safe)**(1/3)

    return rt

def hydrationshell(alpha_t,K3,N3):
    """
    """
    rt = K3 * (1-alpha_t)**N3
    return rt

def f_wc_lothenbach(wc,H,alpha_t):
    """
    functional dependence on w/c in parrot and killoh model
    """
    alpha_cr  = H*wc
    if alpha_t > alpha_cr:
        return (1+3.333*(H*wc-alpha_t))**4
    else:
        return 1

def f_rh(RH):
    """
    functional dependence on relative humidity
    """
    return ((RH-0.55)/0.45)**4

def f_temp(T,T0,Ea):
    """
    functional dependence on temperature
    """
    exp = np.exp
    R = 8.314
    T0  = T + 273.15
    T  =  T + 273.15
    return  exp((-Ea/R) * ((1/T) - (1/T0)))

def overall_rate(t,alpha_t,K1,N1,K2,K3,N3,H,wc,RH,fineness,ref_fineness,T,T0,Ea):
    """
    overall rate of cement clinker phase based on parrot and killoh model
    """
    if alpha_t >= 1: alpha_t = 0.9999
    r1= nuclandgrowth(alpha_t,K1,N1,fineness,ref_fineness)
    r2= diffusion(alpha_t,K2)
    r3= hydrationshell(alpha_t,K3,N3)
    r= min(r1,r2,r3)
    r= r * f_wc_lothenbach(wc,H,alpha_t) * f_rh(RH) * f_temp(T,T0,Ea)
    return r

def parrot_killoh(wc, RH, T, fineness):
    DoH = {}
    for phase in ["C3S","C2S","C3A","C4AF"]:
        ivp_out = solve_ivp(overall_rate, [0,np.max(output_times)], [1e-15],t_eval=output_times
                           ,args=(K1[phase],N1[phase],K2[phase],K3[phase],N3[phase],
                           H[phase],wc,RH,fineness,ref_fineness,T,T0,Ea[phase]))
        DoH[phase] = np.squeeze(ivp_out.y[0])

    return DoH

def plot_bars(results):
    #fix color notation for each phase
    unique_phases = []
    for t in results.keys():
        for phase in results[t].keys():
            if results[t][phase] > 0:
                if phase not in unique_phases:
                    unique_phases.append(phase)
    colors = plt.get_cmap('rainbow')(np.linspace(0., 1, len(unique_phases)+1))
    phase_colors = {}
    for phase,color in zip(unique_phases,colors):
        phase_colors[phase]=color
    i  = 0
    ytks= []
    ytksval = []
    for t in results.keys():
        prev =0
        for key,val in results[t].items():
            if key in unique_phases:
                plt.barh(i, val, left=prev, height=0.5,
                        label=key,color =phase_colors[key])
                prev += val
        ytks.append(i)
        ytksval.append(str(t))
        i+=1
        if t==0:
            plt.legend(ncol=3,bbox_to_anchor=(0, 1.5), loc='upper left')
    plt.yticks(ytks,ytksval)
    plt.xlabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    plt.ylabel("Time [days]")    
    return plt  

# def to_phase_first_dict(results):
#     """
#     convert time first dict to phase first dict for phase diagrams
#     """  
#     for t in results.keys():
#         if t== 0: 
#             out = {}
#             for phase in results[t].keys():
#                 out[phase] = []
#         for phase in results[t].keys():
#             out[phase].append(results[t][phase])
#     for phase in out.keys():
#         out[phase] = np.array(out[phase])
#     return out
def to_phase_first_dict(results):
    out = {}  # Initialize the 'out' dictionary
    for t in results.keys():
        for phase in results[t].keys():
            if phase not in out:  # Initialize the list for this phase if it doesn't exist yet
                out[phase] = []
            out[phase].append(results[t][phase])
    
    # Convert lists to numpy arrays if needed
    for phase in out.keys():
        out[phase] = np.array(out[phase])
    
    return out

# def phase_plot(results,datatype="vfrac"):
#     """
#     plot phase diagrams for visualization
#     """
#     unique_phases=[]
#     for phase in results.keys():
#         if np.count_nonzero(results[phase]) > 0: 
#             unique_phases.append(phase)
#     colors = plt.get_cmap('rainbow')(np.linspace(0., 1, len(unique_phases)+1)) #viridis, rainbow       
#     phase_colors = {}
#     for phase,color in zip(unique_phases,colors):
#         phase_colors[phase]=color
#     prev = 0
#     for phase in unique_phases:
#         plt.fill_between(output_times, results[phase]+prev, prev,color=phase_colors[phase],label=phase)
#         prev += results[phase]
#     plt.legend(ncol=4,bbox_to_anchor=(0, 1.5), loc='upper left')
#     if datatype=="vfrac":
#         plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
#     elif datatype=="masses":
#         plt.ylabel('mass [g/100 g of cement]')

#     plt.xlabel("Time [days]") 
#     plt.xscale("log")
#     if datatype=="vfrac": plt.ylim(0,1.2)
#     plt.xlim((output_times[2]+output_times[1])/2.,output_times[-1])
#     return plt     







# def phase_plot(results, datatype="vfrac"):
#     """
#     Plot phase diagrams for visualization with additional hatches.
#     """
#     unique_phases = []
#     for phase in results.keys():
#         if np.count_nonzero(results[phase]) > 0:
#             unique_phases.append(phase)
    
#     # 扩展的填充图案
#     custom_hatches = ['/', '-', '|', '-', '+', 'x', '|', 'O', '\\', '/', 'H', 'W', 'M', 'N', 'A', 'B', 'C', 'D']
#     colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_phases)))#cividis,tab20
    
#     # 确保填充图案的数量足够多
#     phase_hatches = {phase: custom_hatches[i % len(custom_hatches)] for i, phase in enumerate(unique_phases)}
#     phase_colors = {phase: colors[i] for i, phase in enumerate(unique_phases)}

#     # 提取时间数据的长度
#     num_output_times = len(output_times)

#     prev = 0
#     for phase in unique_phases:
#         plt.fill_between(range(num_output_times), results[phase] + prev, prev, 
#                          color=phase_colors[phase], hatch=phase_hatches[phase], 
#                          edgecolor="black", label=phase, linewidth=0.5)
#         prev += results[phase]

#     # 自定义图例显示填充图案
#     handles = [mpatches.Patch(facecolor=phase_colors[phase], hatch=phase_hatches[phase], 
#                               edgecolor='black', label=phase) for phase in unique_phases]
    
#     plt.legend(
#         handles=handles, 
#         ncol=4, 
#         bbox_to_anchor=(0, 1.5), 
#         loc='upper left', 
#         borderaxespad=0.,
#         fontsize='small',
#         title='Phases',
#         handlelength=3,  # 增加图例图标的长度
#         handletextpad=1.5  # 增加图例项的文本间距
#     )

#     if datatype == "vfrac":
#         plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
#     elif datatype == "masses":
#         # plt.ylabel('Mass [g/100 g of cement]')
#         plt.ylabel('Mass g')
    
#     plt.xlabel("Time [days]") 
    
#     if datatype == "vfrac":
#         plt.ylim(0, 1)
#     plt.xlim(0, 12)  # 设置 x 轴范围
    
#     # 设置 x 轴刻度和标签
#     plt.xticks(ticks=np.arange(13), labels=output_times[:13])
    
#     return plt


# def phase_plot(results, datatype="vfrac"):
#     """
#     Plot phase diagrams for visualization with additional hatches.
#     """
#     unique_phases = []
#     for phase in results.keys():
#         if np.count_nonzero(results[phase]) > 0:
#             unique_phases.append(phase)
#     print(unique_phases)
#     # 扩展的填充图案
#     custom_hatches = ['/', '-', '|', '-', '+', 'x', '|', 'O', '\\', '/', 'H', 'W', 'M', 'N', 'A', 'B', 'C', 'D']
#     colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_phases)))#cividis,tab20
    
#     # 确保填充图案的数量足够多
#     phase_hatches = {phase: custom_hatches[i % len(custom_hatches)] for i, phase in enumerate(unique_phases)}
#     phase_colors = {phase: colors[i] for i, phase in enumerate(unique_phases)}

#     prev = 0
#     for phase in unique_phases:
#         plt.fill_between(output_times, results[phase] + prev, prev, 
#                          color=phase_colors[phase], hatch=phase_hatches[phase], 
#                          edgecolor="black", label=phase, linewidth=0.5)
#         prev += results[phase]

#     # 自定义图例显示填充图案
#     handles = [mpatches.Patch(facecolor=phase_colors[phase], hatch=phase_hatches[phase], 
#                               edgecolor='black', label=phase) for phase in unique_phases]
    
#     plt.legend(
#         handles=handles, 
#         ncol=4, 
#         bbox_to_anchor=(0, 1.5), 
#         loc='upper left', 
#         borderaxespad=0.,
#         fontsize='small',
#         title='Phases',
#         handlelength=3,  # 增加图例图标的长度
#         handletextpad=1.5  # 增加图例项的文本间距
#     )

#     if datatype == "vfrac":
#         plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
#     elif datatype == "masses":
#         # plt.ylabel('Mass [g/100 g of cement]')
#         plt.ylabel('Mass g')
    
#     plt.xlabel("Time [days]") 
    
#     if datatype == "vfrac":
#         plt.ylim(0, 1)
    
#     plt.xscale('log')  # 设置 x 轴为对数刻度
#     # 设置 x 轴刻度和标签
    
#     return plt
    

def phase_plot(results, datatype="vfrac"):
    """
    Plot phase diagrams for visualization with additional hatches.
    """
    unique_phases = []
    
    # 筛选出最大值大于 1 的相
    for phase in results.keys():
        if np.max(results[phase]) > 0.01:
            unique_phases.append(phase)
    
    # print(unique_phases)
    
    # 扩展的填充图案
    custom_hatches = ['/', '-', '|', '-', '+', 'x', '|', 'O', '\\', '/', 'H', 'W', 'M', 'N', 'A', 'B', 'C', 'D']
    colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_phases)))  # cividis, tab20
    
    # 确保填充图案的数量足够多
    phase_hatches = {phase: custom_hatches[i % len(custom_hatches)] for i, phase in enumerate(unique_phases)}
    phase_colors = {phase: colors[i] for i, phase in enumerate(unique_phases)}
    
    prev = 0
    for phase in unique_phases:
        plt.fill_between(output_times, results[phase] + prev, prev, 
                         color=phase_colors[phase], hatch=phase_hatches[phase], 
                         edgecolor="black", label=phase, linewidth=0.5)
        prev += results[phase]

    # 自定义图例显示填充图案
    handles = [mpatches.Patch(facecolor=phase_colors[phase], hatch=phase_hatches[phase], 
                              edgecolor='black', label=phase) for phase in unique_phases]
    
    plt.legend(
        handles=handles, 
        ncol=4, 
        bbox_to_anchor=(0, 1.3), 
        loc='upper left', 
        borderaxespad=0.,
        fontsize='small',
        title='Phases',
        handlelength=3,  # 增加图例图标的长度
        handletextpad=1.5  # 增加图例项的文本间距
    )
    
    if datatype == "vfrac":
        plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    elif datatype == "masses":
        plt.ylabel('Mass g')
    
    plt.xlabel("Time [days]") 
    
    if datatype == "vfrac":
        plt.ylim(0, 1)
    
    # plt.xscale('log')  # 设置 x 轴为对数刻度
    
    # 设置 x 轴刻度和标签
    
    return plt

def co2_phase_plot(results,datatype):
    co2_values=list(results.keys())
    """
    在不同的CO2值下绘制phase的体积分数堆叠图。
    """
    unique_phases = []
    
    # 找出所有不同的相
    if datatype == "vfrac":
        for data in results.values():
            for phase in data.keys():
                if phase not in unique_phases and np.max(data[phase]) > 0.01:
                    unique_phases.append(phase)
    if datatype == "masses":
        for data in results.values():
            for phase in data.keys():
                if phase not in unique_phases and np.max(data[phase]) > 1:
                    unique_phases.append(phase)   
    
    # 自定义填充图案和颜色
    custom_hatches = ['/', '-', '|', '-', '+', 'x', '|', 'O', '\\', '/', 'H', 'W', 'M', 'N', 'A', 'B', 'C', 'D']
    colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_phases)))
    
    # 确保填充图案的数量足够多
    phase_hatches = {phase: custom_hatches[i % len(custom_hatches)] for i, phase in enumerate(unique_phases)}
    phase_colors = {phase: colors[i] for i, phase in enumerate(unique_phases)}
    
    prev = np.zeros(len(co2_values))
    for phase in unique_phases:
        phase_values = [results[co2].get(phase, 0) for co2 in co2_values]
        plt.fill_between(co2_values, prev, prev + np.array(phase_values), color=phase_colors[phase], 
                         hatch=phase_hatches[phase], edgecolor="black", label=phase)
        prev += phase_values
    
    # 自定义图例
    handles = [mpatches.Patch(facecolor=phase_colors[phase], hatch=phase_hatches[phase], 
                              edgecolor='black', label=phase) for phase in unique_phases]
    
    plt.legend(
        handles=handles, 
        ncol=4, 
        bbox_to_anchor=(0, 1.4), 
        loc='upper left', 
        borderaxespad=0.,
        fontsize='small',
        title='Phases',
        handlelength=3,
        handletextpad=1.5
    )

    if datatype == "vfrac":
        plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    elif datatype == "masses":
        plt.ylabel('Mass g')
    
    plt.xlabel("Amount of CO2 [g/100 cement blend]") 
    
    if datatype == "vfrac":
        plt.ylim(0, 1.2)
    plt.show()
    return plt



def run_hydration(clink_phases, wc, CSH2,T, DoH):
   input_file = 'gems_files/CemHyds-dat.lst'
   gemsk = GEMS(input_file)
   all_species = clink_phases.copy()
   all_species["H2O@"] = wc * 100
   # all_species["CO2@"] = recipe['CO2']
   all_species["Gp"]=CSH2
   # print(all_species)
   for name in all_species.keys():
       all_species[name]*=1e-3
   gemsk.T = T + 273.15
   gemsk.add_multiple_species_amt(all_species,units = "kg")
   gemsk.add_species_amt("O2",1e-6) # to reduce stiffness related to redox
   gems_vol_frac = {}
   gems_phase_masses={}
   density = []
   for i in range(len(output_times)):
       if i > 0: gemsk.warm_start()
       for phase in clink_phases:
       #    set lower bound to limit dissolution for unreacted phases
           gemsk.species_lower_bound(phase,clink_phases[phase]*(1-DoH[phase][i])*1e-3,units="kg")
       print("Time-->",str(output_times[i]),":",gemsk.equilibrate())
       gems_phase_masses[output_times[i]]=gemsk.phase_masses
       if output_times[i]==0: init_vol = gemsk.system_volume
       print(init_vol)
       gems_vol_frac[output_times[i]]=gemsk.phase_volumes 
       for key,val in gems_vol_frac[output_times[i]].items():
           gems_vol_frac[output_times[i]][key]/=init_vol
       density.append(gemsk.system_mass/gemsk.system_volume)
   return gems_vol_frac, gems_phase_masses, density

# def run_hydration(clink_phases,remaining_phases, wc, CSH2, T, pk):

#    input_file = 'gems_files/CemHyds-dat.lst'
#    gemsk = GEMS(input_file)
#    # all_species = {**clink_phases, **remaining_phases}
#    all_species = clink_phases.copy()
#    all_species["H2O@"] = wc * 100
#    all_species["Gp"]=CSH2
#    print(all_species)
#    for name in all_species.keys():
#        all_species[name]*=1e-3
#    gemsk.T = T + 273.15
#    gemsk.add_multiple_species_amt(all_species,units = "kg")
#    gemsk.add_species_amt("O2",1e-6) # to reduce stiffness related to redox
#    gems_vol_frac = {}
#    gems_phase_masses={}
#    density = []
#    for i in range(len(output_times)):
#        if i > 0: gemsk.warm_start()
#        for phase in clink_phases:
#        #    set lower bound to limit dissolution for unreacted phases
#            gemsk.species_lower_bound(phase,clink_phases[phase]*(1-DoH[phase][i])*1e-3,units="kg")


#        for component in remaining_phases.keys(): # add SCM
#             if component not in ['clink_phases', 'wc', 'CSH2', 'T', 'RH', 'fineness']:
#                 assert component in SCM.keys()
#                 print(SCM)
#                 try:
#                     value = recipe[component] *1e-3 * alpha(output_times[i], component) *1e-2
#                 except KeyError: # if the component is not in the dic_5pl
#                     value = recipe[component] *1e-3
#                 if SCM[component][1] == 'mol':
#                     gemsk.add_amt_from_formula(SCM[component][0], value, units='mol')
#                 elif SCM[component][1] == 'kg':
#                     gemsk.add_amt_from_formula(SCM[component][0], value, units='kg')
#                 else:
#                     raise UnitsError ('the unit is not taken into account please select "mol" or "kg"')
#        gemsk.supress_phase('gas_gen')
#        equilibration = gemsk.equilibrate()
        
#         #if equilibration fails
#        if equilibration[0] == 'F' or equilibration[0] == 'B': 
#             #print('Equilibration fails for recipe: ', recipe, equilibration, i)
#            return None
#         #if equilibration succeed
#        else: 
#             #suppress the gas phase as it would have enormous volume
#            gemsk.supress_phase('gas_gen')

       
#        print("Time-->",str(output_times[i]),":",gemsk.equilibrate())
#        gems_phase_masses[output_times[i]]=gemsk.phase_masses
#        if output_times[i]==0: init_vol = gemsk.system_volume
#        gems_vol_frac[output_times[i]]=gemsk.phase_volumes 
#        for key,val in gems_vol_frac[output_times[i]].items():
#            gems_vol_frac[output_times[i]][key]/=init_vol
#        density.append(gemsk.system_mass/gemsk.system_volume)
#    return gems_vol_frac, gems_phase_masses, density