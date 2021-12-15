import numpy as np
import os
import sys
import evac_large as ev
#sys.path.insert(1,'/p/project/jias70/jps_jureca/files/jureca_upload')

os.system("mkdir trajectories")
b =  np.array([14.])
start = 0
size = 1
run_jump = 1
pi_round = np.round(np.pi,3)

def ini_list(list):
    list_ini = np.empty(0)
    for element in list:
        list_ini = np.append(list_ini,element)
    return list_ini
    
def list_writer(file,list,list_str):
    if list.shape[0] > 1:
        file.write("{} = np.array([{},".format(list_str,list[0]))
        for element in list[1:-1]:
            file.write("{},".format(element))
        file.write("{}])\n".format(list[-1]))
    else:
        file.write("{} = np.array([{}])\n".format(list_str,list[0]))
    return list
#print(esigma_ini)
#esigma_list = np.array([np.round(np.arange(,0.3,0.05),3)])
#esigma_list = np.round(np.array([[i] for i in np.arange(-1.,1.2,0.2)]),2)
#print(esigma_list)
esigma_list = np.array([[0]])
#d_list = np.array([[0.3],[0.2],[0.1],[0.05]])
d_list = np.array([[0.1]])
#a_list = np.array([[5.],[4.],[3.],[2.5],[2.],[1.5],[1.]])
a_list = np.array([[2.5],[5]])
T_list = np.array([[1.]])

#esigma_ini = np.empty(0)
#for esig in esigma_list:
    #esigma_ini = np.append(esigma_ini,esig)
esigma_ini = ini_list(esigma_list)
a_ini = ini_list(a_list)
d_ini = ini_list(d_list)
T_ini = ini_list(T_list)

print("esigma: ", esigma_ini)
print("d: ", d_ini)
print("a: ", a_ini)
print("T: ", T_ini)


#i_step = 1
#irange = np.arange(i_start,i_final,i_step)
#brange = np.arange(b_min,b_max,b_step)
#print(brange.shape[0])
np.save("ini_b.npy", b)

ev.ini_files(b,start,size,esigma_ini,a_ini,d_ini,T_ini)

mc_i = 1
for T in T_list:
    for a in a_list:
        for d in d_list:
            for esigma in esigma_list:
                for i in np.arange(0,size)[::run_jump]:
                    ranger = int(i + run_jump)
                    if ranger > size:
                        ranger -= ranger - size
                    i_start = str(i)
                    i_final = str(i + run_jump)
                    
                    file = open("mc_l_" + str(mc_i) + ".py", "w")
                    file.write("import sys \n")
                    file.write("sys.path.insert(0,'../') \n")
                    file.write("import evac_large as ev \n")
                    file.write("import numpy as np \n")
                    file.write("i_start = " + i_start + " \n")
                    file.write("i_end = " + i_final + " \n")
                    #file.write("b = np.array([30]) \n")
                    file.write("b = np.array([{}]) \n".format(b[0]))
                    file.write("b = np.array([round(i,3) for i in b]) \n")
                    list_writer(file,esigma,"esigma")
                    list_writer(file,d,"d")
                    list_writer(file,a,"a")
                    list_writer(file,T,"T")


                    
                    """if esigma.shape[0] > 1:
                        file.write("esigma = np.array(["+ str(esigma[0]) + ",")
                        for sigma in esigma[1:-1]:
                            file.write(str(sigma) + ",")
                        file.write(str(esigma[-1]) + "])\n")
                    else:
                        file.write("esigma = np.array(["+str(esigma[0])+"])\n")"""

                    file.write("ev.main(b,i_start,i_end,esigma,a,d,T,i_start) \n")

                    file.close()
                    mc_i += 1
        
        
