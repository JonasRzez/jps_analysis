import sys 
sys.path.insert(0,'../') 
import evac_large as ev 
import numpy as np 
i_start = 0 
i_end = 20 
b = np.array([14.0]) 
b = np.array([round(i,3) for i in b]) 
esigma = np.array([0])
ev.main(b,i_start,i_end,esigma,i_start) 
