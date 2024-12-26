#!/usr/bin/env python3


"""
   This script runs "brownian_beads.py" and records the 
   positions to one file and Rg & Ree to another file.
"""


import brownian_beads as bb
import numpy as np
import os

nbeads_lst = (90, 100)
k_ev = 200
configuration = "circular"
Δt = .0001
t = []
tt = 0
for i in range(bb.N):
    t.append(tt)
    tt += Δt


for nbeads in nbeads_lst:
    print(f"[a] Simulating {configuration} chain: {nbeads} beads, N = {bb.N}, Δt = {Δt}, k_ev = {k_ev}")

    sim = bb.Simulation(nbeads, conf=configuration)
    sim.advance(Δt, κ_ev=k_ev)


    ###################### 'write (wb)' only, no sense to ave pos, disable for multiple...
    print("[b] Saving coordinates")

    array0 = np.array(bb.all_sim_pos)
    with open(f'pos_data/data/circular/k_ev=200/nbeads{nbeads}-k_ev{k_ev}.npy', 'wb') as f:
        np.save(f, array0)
    array1 = np.load(f'pos_data/data/circular/k_ev=200/nbeads{nbeads}-k_ev{k_ev}.npy')

    # make sure arrays are the same
    comparison = (array0 == array1)     # don't need parantheses
    equal_arrays = comparison.all()
    if equal_arrays != True:
        print(equal_arrays)
    ######################


    Ree = []
    for xy in bb.end_to_end:
        Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

    print("[c] Writing Rg & Ree")

    file = open(f"simulated_data/circular/k_ev=200/ave_rg_ree_n{nbeads}-k_ev{k_ev}.dat", "w")  # "a" vs "w"
    for i, s in enumerate(bb.Rg):

        a = str( bb.Rg[i] )
        b = str( Ree[i] )
        c = a + " " + b

        file.write(c + '\n')
    file.close()


    # resetting brownian_bead.py variables for loop
    bb.all_pos_xy = []  
    bb.all_sim_pos = [] 
    bb.end_to_end = []
    bb.Rg = []
    if nbeads != nbeads_lst[-1]:
        print("=========")
#os.system('cmd.exe /c shutdown /s')
