#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

nbeads = 20
k_ev = 200

arr = np.load(f'pos_data/nbeads{nbeads}-k_ev{k_ev}.npy')  
N = len(arr[0])

select_cos_arr = []; 

ı_lst = []
ı = 100; 
while ı < N: 
    temp_x = []
    temp_y = []
    for bead in arr:
        temp_x.append(bead[ı, 0])
        temp_y.append(bead[ı, 1])
    
    vec_arr = np.array([])
    for i in range(1, nbeads):
        vec_arr = np.append(vec_arr, (temp_x[i]-temp_x[i-1], temp_y[i]-temp_y[i-1]))
    vec_arr = np.reshape(vec_arr, (nbeads-1, 2))

    cos_arr = np.array([])
    for i in range(1, nbeads-1):
        cos_γ = np.dot(vec_arr[i],vec_arr[i-1])/(np.linalg.norm(vec_arr[i])*np.linalg.norm(vec_arr[i-1]))
        cos_arr = np.append(cos_arr, np.arccos(cos_γ)*(180/np.pi)) 

    select_cos_arr = np.append(select_cos_arr, cos_arr, axis=0)
    ı_lst.append(ı)
    ı += 100
del ı; 


plt.hist(select_cos_arr, bins=50)
plt.xlabel(r"${\gamma}$", fontsize=20) # \cos{\gamma}
plt.ylabel("frequency", fontsize=20)
plt.tight_layout()
plt.savefig(f"pos_data/pics/k_ev{k_ev}nbeads{nbeads}_histogram.png", bbox_inches='tight')
plt.show()
