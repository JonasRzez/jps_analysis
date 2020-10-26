import numpy as np
import os
import sys
import evac_large as ev
#sys.path.insert(1,'/p/project/jias70/jps_jureca/files/jureca_upload')

os.system("mkdir trajectories")
b =  np.array([50])
start = 0
size = 10
run_jump = 1
esigma = np.array([0.2,0.4,0.5,1.5])
#i_step = 1
#irange = np.arange(i_start,i_final,i_step)
#brange = np.arange(b_min,b_max,b_step)
#print(brange.shape[0])
np.save("ini_b.npy", b)

ev.ini_files(b,start,size,esigma)

mc_i = 1

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
    file.write("b = np.array(["+str(b[0])+"]) \n")
    file.write("b = np.array([round(i,3) for i in b]) \n")
    file.write("esigma = np.array(["+ str(esigma[0]) + ",")
    for sigma in esigma[1:-1]:
    	file.write(str(sigma) + ",")
    file.write(str(esigma[-1]) + "])\n")
    file.write("ev.main(b,i_start,i_end,esigma) \n")

    file.close()
    mc_i += 1
    
    