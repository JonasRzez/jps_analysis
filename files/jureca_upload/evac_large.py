#!/usr/bin/env python
import os
import numpy as np
import itertools
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import random as rand
PATH = os.path.dirname(os.path.abspath("evac_geo_temp.xml"))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(PATH),
    trim_blocks=False)

def csv_writer(file_name,path,df,key,append):
    location = path + file_name
    if os.path.isfile(location) and append:
        load_frame = pd.read_csv(location)
        b_frame = load_frame[key]
        b_np = np.array(b_frame)
        for bi in b_np:
            df = df[df[key] != bi]
        df.to_csv(location, mode = "a", header = False)
    else:
        df.to_csv(location)

def render_template(template_filename, context):
    return TEMPLATE_ENVIRONMENT.get_template(template_filename).render(context)

def Product(variables):
    return list(itertools.product(*variables))

def create_inifile(geo_name, geo_name_ini, cross_new, ini_list, ini_folder_list , location, folder,  stepsize, fps, l0_list, t_max,periodic, rand_mot, mot_hm_lm,ini_i,model,dr):
    
   #location = [path + "evac_traj_" + str(2 * bi)[0] + str(2*bi)[-1] + ".txt" for bi in b_list]
    #print(geo_name,location)
    
    for var,fname,l0 in zip(cross_new, geo_name,l0_list):
        context= {'b': var[1],'l':l0, 'll': 1.5 * l0,"wedge":var[1] - 0.4}
        fname = folder + fname
        if os.path.isfile(fname) == False:
            #print("output geo: ",  fname)
            print(os.system("pwd"))
            with open(fname,'w') as f:
                xml = render_template('evac_geo_temp.xml', context)
                f.write(xml)
    #print("shapes " , np.array(geo_name_ini).shape,np.array(location).shape,np.array(ini_list).shape,np.array(cross_new).shape)
    
    for geo,loc,fname_ini,var,l0,ini_folder_i in zip(geo_name_ini,location,ini_list,cross_new,l0_list,ini_folder_list):
        print("ini_folder_i = ", ini_folder_i)
        print("fname iini = " , fname_ini)
        print("folder = ", folder)
        seed = int(rand.uniform(0, 1485583981))
        stepsize = modelStepSize(var[7])
        if rand_mot:
            ini_context = {'geo':geo,'location':loc,'b':var[1],'l':l0, 'll':  1.5 * l0,'seed':seed,'stepsize':stepsize,'fps':fps,'N_ped_lm':int(var[5] * (1 - var[6])),'N_ped_hm': int(var[5] * var[6]),'esig':var[0],'t_max':t_max,'periodic':periodic,'v0':var[2],'T_lm':mot_hm_lm[1],'T_hm':mot_hm_lm[0],'rsigma':var[8],'r':var[9],'a':var[10],'d':var[11], "aviod_wall":var[12],'output': "../../files/jureca_upload/" + folder + ini_folder_i }
        else:
            ini_context = {'geo':geo,'location':loc,'b':var[1],'l':l0, 'll':  1.5 * l0,'seed':seed,'stepsize':stepsize,'fps':fps,'N_ped_lm':int(var[5] * (1 - var[6])),'N_ped_hm': int(var[5] * var[6]),'esig':var[0],'t_max':t_max,'periodic':periodic,'v0':var[2],'T_lm':var[3],'T_hm':var[3],'rsigma':var[8],'r':var[9],'r2':var[9] + dr,'a':var[10],'d':var[11], "aviod_wall":var[12],'output': "../../files/jureca_upload/" + folder + ini_folder_i}

        #l = 0
        #l_add = round(l0/6,2)
        """for i in range(1,7): #devides the room into 6 slices to distribute pedestrians
            l += l_add
            ini_context["l" + str(i)] = l"""
        #print("ini_context = ", ini_context)
        fname_ini = folder + fname_ini
        #print("output ini: ", fname_ini)
        with open(fname_ini, 'w') as f:
            xml_ini = render_template('evac_ini_temp_diff_{}.xml'.format(modelStr(var[7])), ini_context)
            f.write(xml_ini)
            
def modelStepSize(model_number):
    if model_number == 0:
        return 0.05
    if model_number == 1:
        return 0.005
    print("ERROR: model number does not exist")
    
            
def b_data_name(b,dig):
    b = round(b,dig)
    str_b = ''
    for let in str(b):
        if (let in '.'):
            str_b += '_'
        else:
            str_b += let
    return str_b
    
def test_var_fetch(test_var):
    if test_var == "d":
        var_i = 11
    if test_var == "a":
        var_i = 10
    if test_var == "r_i":
        var_i = 9
    if test_var == "rsigma":
        var_i = 8
    if test_var == "model":
        var_i = 7
    if test_var == "mot_frac":
        var_i = 6
    if test_var == "N_ped":
        var_i = 5
    if test_var == "rho":
        var_i = 4
    if test_var == "T":
        var_i = 3
    if test_var == "v0":
        var_i = 2
    if test_var == "b":
        var_i = 1
    if test_var == "esigma":
        var_i = 0

    else:
        print("WARNING: test var not included")
    return var_i

def ini_bool():
    sec_test_var = True
    append = False
    rho_ini_rand = False
    rand_mot = False
    polydispers = True
    return sec_test_var,append,rho_ini_rand,rand_mot,polydispers
def test_var_ini():
    dig = 3  # for the rounding of digits
    test_var = "T"
    test_var2 = "d"
    motivation = "_lm_"
    var_i = test_var_fetch(test_var)

    return dig,test_var,test_var2,var_i,motivation

def var_ini(i_start,i_end,esigma,a_array,d_array,T):
    #rho_min = 2.0
    #rho_max = 3.0
    rho_min = 4.
    rho_max = 4.
    #rho_ini = np.array([(rho_min + rho_max) / 2])
    rho_ini = np.array([4.3])

    #rsigma = np.array([0.0,0.01,0.02,0.04,0.05,0.06,0.07,0.08,0.09,0.1])
    #rsigma = np.array([0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.09])
    rsigma = np.array([0.0])
    #rsigma = np.array([0.0])
    #r_array = [0.15,0.16,0.17,0.18,0.19,0.2]
    r_array = np.array([0.16])
    #a_array = np.array([5,2.5,1.5,1])
    #a_array = np.array([2.5])
    #d_array = np.array([0.1,0.05])
    #d_array = np.array([0.2,0.15,0.1,0.05,0.01])
    #r_array = np.array([0.16,0.165,0.17,0.175,0.18])
    #r_array = np.array([0.17])
    #rsigma = np.array([0.0])
    #rsigma = np.array([0.5,0.8,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,4.0])
    
    #T = np.array([0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    #T = np.array([0.1,0.15,0.2,0.25,0.3])
    #T = np.array([0.25])
    
    #T = np.array([0.6])
    #mot_frac = np.array([1.,0.8,0.6,0.5,0.4])
    mot_frac = np.array([1.])

    #poly = np.array([1.0,0.8,0.4,0.2,0.0])
    avoid_wall = np.array([0.0,0.0])
    v0 = np.array([1.34])
    #esigma = np.array([0.7])
    N_ped = np.array([20])
    model = [0]
    #r = 0.16
    #r_array = np.array([0.16])
    dr = 0.0

    fps = 2
    stepsize = 0.005
    #N_ped = 55
    N_runs = 1
    t_max = 10
    
    periodic = 1
    return rho_ini, T, v0, esigma, fps, stepsize, N_ped, i_start, i_end, t_max, periodic, N_runs, rho_min, rho_max, avoid_wall, mot_frac, model, rsigma, r_array,dr, a_array, d_array

def ini_cross(cross_variable,shape_var,rho_min,rho_max,rho_ini_rand):
    cross_new = np.empty([shape_var,
                          13])  # length has to be the number of variables given currently 0:esigma,1:b, 2:v0, 3:T , 4:rho_ini: 5 N_ped 6: mot_frac, 7: model, 8:rsigma 9:r, 10:a, 11:d
    avoid_wall = np.array([0.4])
    #T = np.array([1.0, 0.9,0.8,0.7])
    #print(cross_new)
    for i in range(shape_var):
        cross_i = cross_variable[i]
        cross_i = np.append(cross_i,avoid_wall[0])
        print(cross_i.shape, cross_i)
        cross_new[i] = cross_i  # vraiable are 0:esigma, 1:b, 2:v0, 3:T, 4:rho_ini, 5:N_ped, 6:mot_frac, 7:model, 8:rsigma, 9:r ,10:a,11:d, 12:wall avoidance
    return cross_new
    
def modelStr(model):
    model_str = " "
    model_list = [0,1]
    if int(model) == 0:
        model_str = "velocity"
    if int(model) == 1:
        model_str = "force"
    elif model not in model_list:
        print("ERROR: Model not in list of Models. Model ID {}  does not exist".format(int(model)))
    return model_str

def ini_csv(path,N_runs,fps,N_ped,t_max,periodic,var_i,test_var,sec_test_var,test_var2,b,append,v0,T,esigma,rho_ini,ini_folder_list,data_folder,cross_new,model,rsigma,r_array,a_array,d_array,variables):
    # <saving all variable information into csv files>

    variables_df = pd.DataFrame(
        {'r_i': [r_array], 'fps': [fps], 'N_runs': [N_runs], 'N_ped': [N_ped[0]], 't_max': [t_max], 'periodic': [periodic],
         'test_var': [var_i], 'test_str': [test_var]})
    if sec_test_var:
        var_j = test_var_fetch(test_var2)
        variables_df['test_var2'] = var_j
        variables_df['test_str2'] = test_var2
    variables_df.to_csv(path + "variables_list.csv", mode="w")
    b_df = pd.DataFrame({'b': b})
    csv_writer("b_list.csv", path, b_df, 'b', append)

    v_df = pd.DataFrame({'v0': v0})
    csv_writer("v0_list.csv", path, v_df, 'v0', append)

    T_df = pd.DataFrame({'T': T})
    csv_writer("T_list.csv", path, T_df, 'T', append)

    N_df = pd.DataFrame({'N_ped': N_ped})
    csv_writer("N_ped_list.csv", path, N_df, 'N_ped', append)

    esig_df = pd.DataFrame({'esigma': esigma})
    csv_writer("esigma_list.csv", path, esig_df, 'esigma', append)

    rho_df = pd.DataFrame({'rho_ini': rho_ini})
    csv_writer("rho_list.csv", path, rho_df, 'rho_ini', append)

    mot_frac_df = pd.DataFrame({'rho_ini': rho_ini})
    csv_writer("mot_frac_list.csv", path, mot_frac_df, 'mot_frac', append)
    
    model_df = pd.DataFrame({'model': model})
    csv_writer("model_list.csv", path, model_df, 'model', append)
    
    rsigma_df = pd.DataFrame({'rsigma': rsigma})
    csv_writer("rsigma_list.csv", path, rsigma_df, 'rsigma', append)
    
    r_df = pd.DataFrame({'r_i': r_array})
    csv_writer("r_i_list.csv", path, r_df, 'r_i', append)
    
    a_df = pd.DataFrame({'a': a_array})
    csv_writer("a_list.csv", path, a_df, 'a', append)
    
    d_df = pd.DataFrame({'d': d_array})
    csv_writer("d_list.csv", path, a_df, 'd', append)

    folder_path_list = np.array([data_folder + "/" + ini for ini in ini_folder_list])

    T_list = np.array([var[3] for var in cross_new])
    esig_list = np.array([var[0] for var in cross_new])
    v0_list = np.array([var[2] for var in cross_new])
    b_list = np.array([var[1] for var in cross_new])
    rho_list = np.array([var[4] for var in cross_new])
    N_ped_list = np.array([var[5] for var in cross_new])
    mot_frac_list = np.array([var[6] for var in cross_new])
    model_list = np.array([var[7] for var in cross_new])
    rsigma_list = np.array([var[8] for var in cross_new])
    r_list = np.array([var[9] for var in cross_new])
    a_list = np.array([var[10] for var in cross_new])
    d_list = np.array([var[11] for var in cross_new])


    #print("shapes of list = ", T_list.shape,esig_list.shape,v0_list.shape,rho_list.shape,b_list.shape,np.array(ini_folder_list).shape)
    folder_df = pd.DataFrame(
        {'ini_folder': ini_folder_list, 'b': b_list, 'v0': v0_list, 'T': T_list, 'rho': rho_list, 'esigma': esig_list,'N_ped':N_ped_list, 'mot_frac':mot_frac_list, 'model':model_list,'rsigma':rsigma_list,'r_i':r_list,'a':a_list,'d':d_list})

    # folder_df.to_csv( path + "folder_list.csv")
    csv_writer("folder_list.csv", path, folder_df, test_var, append)

    path_df = pd.DataFrame({'path': [path]})
    path_df.to_csv("path.csv")

    np.save(path + "cross_var.npy", cross_new)
    np.save(path + "var.npy", variables)


# </saving information>

def ini_files(b,i_start,i_end,esigma,a_array,d_array,T):
    sec_test_var, append, rho_ini_rand, rand_mot,polydispers = ini_bool()
    dig, test_var, test_var2, var_i, motivation = test_var_ini()
    rho_ini, T, v0, esigma, fps, stepsize, N_ped, i_start, i_end, t_max, periodic, N_runs, rho_min, rho_max, wall_avoidance, mot_frac, model, rsigma,r_array,dr,a_array,d_array = var_ini(i_start,i_end,esigma,a_array,d_array,T)
    b = b/2

    variables = np.array([esigma, b, v0, T, rho_ini,N_ped, mot_frac,model,rsigma,r_array,a_array,d_array])
    cross_variable = np.array(list(itertools.product(*variables)))
    shape_var = cross_variable.shape[0]

    cross_new = ini_cross(cross_variable, shape_var, rho_min, rho_max, rho_ini_rand)

    path, data_folder, traj_folder = ini_traj_folder(motivation, N_ped, t_max,  fps, test_var)

    ini_folder_list = ini_folder_fetch(cross_new, dig, motivation, N_ped, periodic, traj_folder, data_folder, t_max)

    ini_csv(path, N_runs, fps, N_ped, t_max, periodic, var_i, test_var, sec_test_var, test_var2, b, append, v0, T,
            esigma, rho_ini, ini_folder_list, data_folder, cross_new,model,rsigma,r_array,a_array,d_array, variables)

def ini_folder_fetch(cross_new,dig,motivation,N_ped,periodic,traj_folder,data_folder,t_max):
    ini_folder_list = []
    for var in cross_new:
        ini_folder_name = "ini_" + b_data_name(var[1], dig) + motivation + str(N_ped[0]) + "_esigma_" + b_data_name(var[0],
                                                                                                                 dig) + "_tmax_" + str(
            t_max) + "_periodic_" + str(periodic) + "_v0_" + b_data_name(var[2], dig) + "_T_" + b_data_name(var[3],
                                                                                                            dig) + "_rho_ini_" + b_data_name(
            var[4], dig) + "_Nped_" + b_data_name(var[5],dig) + "_motfrac_" + b_data_name(var[6],dig) + "_model_" + modelStr(var[7]) + "_rsigma_" + b_data_name(var[8],dig) + "_r_" + b_data_name(var[9],dig) + "_a_" + b_data_name(var[10],dig) + "_d_" + b_data_name(var[11],dig)
        mkdir_path = traj_folder + data_folder + "/" + ini_folder_name
        print("os path is file = ", os.path.isdir(mkdir_path))

        if 1 - os.path.isdir(mkdir_path):
            os.system("mkdir " + mkdir_path)
        ini_folder_list.append(ini_folder_name)
    return ini_folder_list

def ini_traj_folder(motivation,N_ped,t_max,fps,test_var):
    traj_folder = "trajectories/"
    data_folder = "ini" + motivation + "N_ped" + str(N_ped[0]) + "_tmax" + str(t_max) +  "_fps_" + str(fps) + "_testvar_" + test_var
    print_level = "--log-level [debug, info, warning, error, off]"
    if 1 - os.path.isfile(traj_folder + data_folder):
        print(traj_folder + data_folder)
        os.system("mkdir " + traj_folder + data_folder)
    path = traj_folder + data_folder + "/"

    return path,data_folder,traj_folder

def main(b,i_start,i_end,esigma,a_array,d_array,T,ini_i):
    sec_test_var, append, rho_ini_rand,rand_mot,polydispers = ini_bool()
    dig,test_var,test_var2,var_i,motivation = test_var_ini()
    rho_ini,T , v0,esigma, fps, stepsize, N_ped, i_start, i_end, t_max, periodic, N_runs,rho_min,rho_max,avoid_wall,mot_frac, model, rsigma, r_array, dr, a_array,d_array = var_ini(i_start,i_end,esigma,a_array,d_array,T)
    run_total = N_runs * b.shape[0]
    b = b / 2
    mot_hm_lm = [0.1,1.3]

    variables = np.array([esigma,b,v0,T,rho_ini,N_ped,mot_frac,model,rsigma,r_array,a_array,d_array])
    #print("variables = " ,variables)
    cross_variable = np.array(list(itertools.product(*variables)))
    shape_var = cross_variable.shape[0]
    print("shape_var = ", shape_var)
    #print("corr variable = ", cross_variable)
    cross_new = ini_cross(cross_variable, shape_var, rho_min, rho_max, rho_ini_rand)

    path, data_folder, traj_folder = ini_traj_folder(motivation,N_ped,t_max,fps,test_var)

    ini_folder_list = ini_folder_fetch(cross_new, dig, motivation, N_ped, periodic, traj_folder, data_folder, t_max)
    #ini_folder_list = pd.read_csv(path + "folder_list.csv")["ini_folder"].values
    geo_name = [data_folder + "/" + ini_folder + "/" + "geo_" + b_data_name(var[1],dig) + ".xml" for var, ini_folder in zip(cross_new, ini_folder_list)]
    geo_name_ini = ["geo_" + b_data_name(var[1],dig) + ".xml" for var in cross_new]

    ini_list = [data_folder + "/" + ini_folder + "/" + "ini_" + b_data_name(var[1],dig)+ str(ini_i) + ".xml" for var, ini_folder in zip(cross_new, ini_folder_list)]
    output_list = [data_folder + "/" + ini_folder for var, ini_folder in zip(cross_new, ini_folder_list)]
    for i in range(i_start,i_end):
        print("iteration = " , i)
        location = ["evac_traj_" + b_data_name(2 * var[1],dig) + "_" + str(i) + ".txt" for var, folder_name in zip(cross_new,ini_folder_list)]
        #print("location = ", location)
        #rho_ini_rand = [rand.uniform(3,4) for i in cross_new]
        if rho_ini_rand:
            for i in range(shape_var):
                cross_new[i][4] = rand.uniform(rho_min, rho_max)

        l0_list = [round(var[5]/(2 * var[4]*var[1]),2) for var in cross_new]
        #l0_list = [l0 if l0 > 7. else 7. for l0 in l0_list]
        create_inifile(geo_name,geo_name_ini,cross_new,ini_list ,output_list,location,traj_folder,stepsize,fps,l0_list,t_max,periodic,rand_mot,mot_hm_lm,ini_i,model,dr)
        os.system("pwd")
        os.chdir("../../build/bin")
        run_count = 0
        
        for ini, ini_folder in zip(ini_list,ini_folder_list):
            jps = "./jpscore ../../files/jureca_upload/" + traj_folder + ini
            #print(jps)
            os.system(jps)
            run_count += 1
        os.chdir("../../files/jureca_upload")

if __name__ == "__main__":
    main(b,i_start,i_end,esigma,a_array,d_array,T,ini_i)


