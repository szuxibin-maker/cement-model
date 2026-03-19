import pandas as pd
import numpy as np
from scipy.stats import qmc
from scipy.stats.qmc import LatinHypercube
from run.hydration import to_phase_first_dict
from util.final_hydration import final_hydration

def convert_nparray_to_string_list(nparray : np.ndarray):
    """ Convert an np.array into a list of strings with same numerical values"""
    list_of_str = []
    for i in nparray:
        list_of_str.append(str(i))
    return list_of_str

def add_new_material(cook_book : dict, material : str, default_value : float, bounds = False):
    """ Add new material to the cookbook with bounds """
    dic = {}
    dic[material] = default_value
    if bounds != False:
        dic['bounds'] = bounds
    cook_book[material] = dic

def create_grid(bound : tuple):
    """ function to create a grid of values """
    n_sample = 28
    start = bound[0]
    end = bound[1]
    grid = np.linspace(start,end, num = n_sample)
    return grid

def create_grid_from_cook_book(cook_book : dict):
    """ Create a grid in a dictionary form from a cook_book"""
    grid = {}
    for param in cook_book.keys():
        if param == 'clink_phase':
            for key in cook_book[param].keys():
                if key != 'bounds':
                    grid[key] = create_grid(cook_book[param]['bounds'][key]) 
        else:
            if 'bounds' in cook_book[param].keys():
                grid[param] = create_grid(cook_book[param]['bounds'])
            else:
                grid[param] = cook_book[param][param]
    return grid

def get_bounds_from_cook_book(cook_book : dict):
    """ Get the bounds from cook_book and the sobol sequence dimension """
    bound = {}
    count_dim = 0 #to compute the dimension of the sobol sequence
    for param in cook_book.keys():
        if param == 'clink_phase':
            for key in cook_book[param].keys():
                if key != 'bounds':
                    bound[key] = cook_book[param]['bounds'][key]
                    count_dim += 1
        else:
            if 'bounds' in cook_book[param].keys():
                bound[param] = cook_book[param]['bounds']
                count_dim += 1
    return bound, count_dim

def normalize_sobol(sobol_seq : np.ndarray, bound: tuple):
    """ Rescale the sobol sequence to the wanted interval caracterised by bound """
    return (bound[1]-bound[0])*sobol_seq + bound[0]

def create_sobol_sequence(cook_book : dict, n_sample: int):
    """ Create a sobol sequence of length 2^n_sample from the cook_book """
    sobol = {}
    bound, sobol_dim = get_bounds_from_cook_book(cook_book)
    assert sobol_dim < 120
    sampler = qmc.Sobol(sobol_dim, seed = 42)
    sample = sampler.random_base2(n_sample)
    count = 0
    for param in cook_book.keys():
        if param == 'clink_phase':
            for key in cook_book[param].keys():
                if key != 'bounds':
                    sobol[key] = normalize_sobol(sample[:,count], bound[key])
                    count += 1
        else:
            if 'bounds' in cook_book[param].keys():
                sobol[param] = normalize_sobol(sample[:,count], bound[param])
            else:
                sobol[param] = cook_book[param][param]
    return sobol

def normalize_samples(samples, bounds):
    """ Normalize samples according to the specified bounds """
    return bounds[0] + samples * (bounds[1] - bounds[0])


def create_lhs_samples(cook_book: dict, n_sample: int):
    """ Create Latin hypercube samples from the cook_book """
    lhs_samples = {}
    bound, lhs_dim = get_bounds_from_cook_book(cook_book)
    assert lhs_dim < 120

    # Generate Latin hypercube samples for each parameter
    for param in cook_book.keys():
        if param == 'clink_phase':
            for key in cook_book[param].keys():
                if key != 'bounds':
                    lhs_sampler = LatinHypercube(d=1)
                    lhs_samples[key] = normalize_samples(lhs_sampler.random(n_sample), bound[key])
        else:
            if 'bounds' in cook_book[param].keys():
                lhs_sampler = LatinHypercube(d=1)
                lhs_samples[param] = normalize_samples(lhs_sampler.random(n_sample), bound[param])
            else:
                lhs_samples[param] = cook_book[param][param]

    return lhs_samples


# def create_recipe(data : pd.DataFrame):
#     """ Create recipes that the hydration model understand from the data """
#     recipe = {}
#     clinker_phases = {}
#     sum_clinker_phases = data.loc[['C3S', 'C2S', 'C3A', 'C4AF', 'CSH2']].sum()
    
#     clinker_phases["C3S"]=data.loc['C3S'] *100 / sum_clinker_phases #alite
#     # divide by the square of the sum to make the sum of all constituants = 100 as in cemgems
#     clinker_phases["C2S"]=data.loc['C2S'] *100 / sum_clinker_phases #belite
#     clinker_phases["C3A"]=data.loc['C3A'] *100 / sum_clinker_phases #tricalcium aluminate
#     clinker_phases["C4AF"]=data.loc['C4AF'] *100 / sum_clinker_phases #calcium aluminoferrite
#     recipe['clink_phases'] = clinker_phases
#     recipe['CSH2'] = data.loc['CSH2']*100 / sum_clinker_phases
#     assert (99 <clinker_phases["C3S"] + clinker_phases["C2S"] + clinker_phases["C3A"] + clinker_phases["C4AF"] + recipe["CSH2"] < 101)
#     sum_SCM = data.loc[~data.index.isin(['C3S','C2S','C3A','C4AF', 'CSH2','wc', 'T', 'RH', 'fineness'])].sum()
#     for add in data.index:
#         if add not in ['C3S','C2S','C3A','C4AF', 'CSH2']:
#             recipe[add] = data.loc[add]
#     recipe['wc'] = data.loc['wc']*(sum_clinker_phases+sum_SCM)
#     recipe['wc_ratio'] = data.loc['wc']
#     return recipe

def create_recipe(data : pd.DataFrame):
    """ Create recipes that the hydration model understand from the data """
    recipe = {}
    clinker_phases = {}
    sum_clinker_phases = data.loc[['C3S', 'C2S', 'C3A', 'C4AF', 'CSH2']].sum()
    
    clinker_phases["C3S"]=data.loc['C3S']  #alite
    # divide by the square of the sum to make the sum of all constituants = 100 as in cemgems
    clinker_phases["C2S"]=data.loc['C2S'] #belite
    clinker_phases["C3A"]=data.loc['C3A']#tricalcium aluminate
    clinker_phases["C4AF"]=data.loc['C4AF']  #calcium aluminoferrite
    recipe['clink_phases'] = clinker_phases
    recipe['CSH2'] = data.loc['CSH2']
    recipe["CO2"] = data.loc['CO2']
    sum_SCM = data.loc[~data.index.isin(['C3S','C2S','C3A','C4AF', 'CSH2','wc', 'T', 'RH', 'fineness','CO2'])].sum()
    for add in data.index:
        if add not in ['C3S','C2S','C3A','C4AF', 'CSH2','CO2']:
            recipe[add] = data.loc[add]
    recipe['wc'] = data.loc['wc']*100
    # recipe['wc_ratio'] = data.loc['wc']
    return recipe


def complete_hydration_from_data(data : pd.DataFrame, output_materials : list, print_error : bool = False):
    """ Compute the hydration from the data generated and compile it with materials of interest contained in output_materials     list """
    vol_frac_dic = {}
    vol_frac_dic2 = {}
    index_to_drop = []
    gems_phase_amounts_dic = {}
    gems_phase_amounts_dic2 = {}
    index_to_drop_gems_phase_amounts = []
    gems_B_dic = {}
    gems_pE_dic = {}
    gems_pH_dic = {}
    gems_density_dic = {}
    gems_aq_volume_frac_dic = {}
    

    
    for i in range(data.shape[0]):
        recipe = create_recipe(data.iloc[i])
        #print(recipe)
        gems_vol_frac , gems_phase_amounts, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac,element_CSHQ,element_aq,other_vol_delete_gas = final_hydration(recipe)
       
        # gems_vol_frac= final_hydration(recipe)
        # print(gems_B)
        # print(gems_phase_amounts)
        # print(gems_vol_frac)
        if gems_vol_frac == None:
            index_to_drop.append(i)
            if print_error == True:
                pass
                #print('Equilibration failed for recipe ', recipe) 
                #print('----------- New Try ------------')
                #gems_vol_frac = final_hydration(recipe)
                #if gems_vol_frac == None:
                #    print('-------- No Equilibration ---------')
        else:
            pd_vol_frac = pd.Series(gems_vol_frac)
            other_material = pd_vol_frac.loc[~pd_vol_frac.index.isin(output_materials)].sum() 
            pd_vol_frac = pd_vol_frac.loc[output_materials] #material of interest
            pd_vol_frac['other_material'] = other_material
            vol_frac_dic[str(i)] = pd_vol_frac
            vol_frac_dic2[str(i)] = pd.Series(gems_vol_frac)
            
        # if gems_phase_amounts == None:
        #     index_to_drop_gems_phase_amounts.append(i)
        #     if print_error == True:
        #         pass
                    
        # else:
            pd_gems_phase_amounts = pd.Series(gems_phase_amounts)/ 0.001
           #  print('pd_gems_phase_amounts',pd_gems_phase_amounts)
           # print('***start')
           # print("end***")
            other_material = pd_gems_phase_amounts.loc[~pd_gems_phase_amounts.index.isin(output_materials)].sum() #sum of all materials that are not of interest
            
            pd_gems_phase_amounts = pd_gems_phase_amounts.loc[output_materials] #material of interest
            pd_gems_phase_amounts['other_material'] = other_material
            gems_phase_amounts_dic[str(i)] = pd_gems_phase_amounts
            gems_phase_amounts_dic2[str(i)] = pd.Series(gems_phase_amounts)
            gems_phase_amounts = gems_phase_amounts_dic
            
            pd_gems_B = pd.Series(gems_B)
            gems_B_dic[str(i)] = pd_gems_B
            gems_B = gems_B_dic
            
            pd_gems_pE = pd.Series(gems_pE)
            gems_pE_dic[str(i)] = pd_gems_pE
            gems_pE = gems_pE_dic

            pd_gems_pH = pd.Series(gems_pH)
            gems_pH_dic[str(i)] = pd_gems_pH
            gems_pH = gems_pH_dic
            
            pd_gems_density = pd.Series(gems_density)
            gems_density_dic[str(i)] = pd_gems_density
            gems_density = gems_density_dic
            
            pd_gems_aq_volume_frac = pd.Series(gems_aq_volume_frac)
            gems_aq_volume_frac_dic[str(i)] = pd_gems_aq_volume_frac
            gems_aq_volume_frac = gems_aq_volume_frac_dic
    
    return index_to_drop, pd.DataFrame(vol_frac_dic, index = pd_vol_frac.keys()), vol_frac_dic2, gems_phase_amounts, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac


def complete_hydration_CO2(data : pd.DataFrame,  print_error : bool = False):
    """ Compute the hydration from the data generated and compile it with materials of interest contained in output_materials     list """
    vol_frac_dic_filter = {}
    vol_frac_dic_all = {}
    index_to_drop = []
    gems_phase_amounts_dic_filter = {}
    gems_phase_amounts_dic_all = {}
    index_to_drop_gems_phase_amounts = []
    gems_B_dic = {}
    gems_pE_dic = {}
    gems_pH_dic = {}
    gems_density_dic = {}
    gems_aq_volume_frac_dic = {}
    gems_other_vol_delete_gas_dic = {}
    element_CSHQ_dic = {}
    
    for i in range(data.shape[0]):
        recipe = create_recipe(data.iloc[i])
        #print(recipe)
        gems_vol_frac , gems_phase_amounts, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac,element_CSHQ,element_aq,other_vol_delete_gas = final_hydration(recipe)
       
        # gems_vol_frac= final_hydration(recipe)
        # print(gems_B)
        # print(gems_phase_amounts)
        # print(gems_vol_frac)
        if gems_vol_frac == None:
            index_to_drop.append(i)
            if print_error == True:
                pass
                print('Equilibration failed for recipe ', recipe) 
                print('----------- New Try ------------')
                gems_vol_frac = final_hydration(recipe)
                if gems_vol_frac == None:
                   print('-------- No Equilibration ---------')
        else:
            pd_vol_frac = pd.Series(gems_vol_frac)
            # other_material = pd_vol_frac.loc[~pd_vol_frac.index.isin(output_materials)].sum() 
            # pd_vol_frac = pd_vol_frac.loc[output_materials] #material of interest
            # pd_vol_frac['other_material'] = other_material
            # vol_frac_dic_filter[str(i)] = pd_vol_frac
            vol_frac_dic_all[str(i)] = pd.Series(gems_vol_frac)
            
        if gems_phase_amounts == None:
            index_to_drop_gems_phase_amounts.append(i)
            if print_error == True:
                pass
                    
        else:
            # pd_gems_phase_amounts = pd.Series(gems_phase_amounts)/ 0.001
           #  print('pd_gems_phase_amounts',pd_gems_phase_amounts)
           # print('***start')
           # print("end***")
            # other_material = pd_gems_phase_amounts.loc[~pd_gems_phase_amounts.index.isin(output_materials)].sum() #sum of all materials that are not of interest
            
            # pd_gems_phase_amounts = pd_gems_phase_amounts.loc[output_materials] #material of interest
            # pd_gems_phase_amounts['other_material'] = other_material
            # gems_phase_amounts_dic_filter[str(i)] = pd_gems_phase_amounts
            gems_phase_amounts_dic_all[str(i)] = pd.Series(gems_phase_amounts)/ 0.001
            # gems_phase_amounts = gems_phase_amounts_dic
            
            pd_gems_B = pd.Series(gems_B)
            gems_B_dic[str(i)] = pd_gems_B
            gems_B = gems_B_dic
            
            pd_gems_pE = pd.Series(gems_pE)
            gems_pE_dic[str(i)] = pd_gems_pE
            gems_pE = gems_pE_dic

            pd_gems_pH = pd.Series(gems_pH)
            gems_pH_dic[str(i)] = pd_gems_pH
            gems_pH = gems_pH_dic
            
            pd_gems_density = pd.Series(gems_density)
            gems_density_dic[str(i)] = pd_gems_density
            gems_density = gems_density_dic
            
            pd_gems_aq_volume_frac = pd.Series(gems_aq_volume_frac)
            gems_aq_volume_frac_dic[str(i)] = pd_gems_aq_volume_frac
            gems_aq_volume_frac = gems_aq_volume_frac_dic

            pd_gems_other_vol_delete_gas = pd.Series(other_vol_delete_gas)
            gems_other_vol_delete_gas_dic[str(i)] = pd_gems_other_vol_delete_gas
            gems_other_vol_delete_gas_volume = gems_other_vol_delete_gas_dic
            

            element_CSHQ = pd.Series(element_CSHQ)
            element_CSHQ_dic[str(i)] = element_CSHQ
            element_CSHQ = element_CSHQ_dic
    
    # return index_to_drop, pd.DataFrame(vol_frac_dic_filter, index = pd_vol_frac.keys()),vol_frac_dic_all, gems_phase_amounts_dic_filter,  gems_phase_amounts_dic_all, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac,element_CSHQ,element_aq
    
    return index_to_drop, pd.DataFrame(vol_frac_dic_all, index = pd_vol_frac.keys()),vol_frac_dic_all,  gems_phase_amounts_dic_all, gems_B,gems_pE,gems_pH,gems_density,gems_aq_volume_frac,element_CSHQ,element_aq,gems_other_vol_delete_gas_volume