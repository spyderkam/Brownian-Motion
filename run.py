#!/usr/bin/env python3


"""
   This script runs "brownian_bead" and records the positions
   to one file and Rg & Ree to another file.
"""


from brownian_bead import *

nbeads = 20
k_ev = 200

Δt = .0001
t = []
tt = 0
for i in range(N):
    t.append(tt)
    tt += Δt

print(f"[a] Simulating ideal chain: {nbeads} beads, N = {N}, Δt = {Δt}, k_ev = {k_ev}")

sim = Simulation(nbeads)
sim.advance(Δt, κ_ev=k_ev)


########################### 'write' only, no sense to ave pos, disable for multiple...
print("[b] Saving coordinates")

array0 = np.array(all_sim_pos)
with open(f'pos_data/nbeads{nbeads}-k_ev{k_ev}.npy', 'wb') as f:
    np.save(f, array0)
array1 = np.load(f'pos_data/nbeads{nbeads}-k_ev{k_ev}.npy')


# make sure arrays are the same
comparison = (array0 == array1)     # don't need parantheses
equal_arrays = comparison.all()
if equal_arrays != True:
    print(equal_arrays)
###########################


Ree = []
for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

print("[c] Writing Rg & Ree")

file = open(f"simulated_data/ave_rg_ree_n{nbeads}-k_ev{k_ev}.dat", "w")  # "a" vs "w"
for i, s in enumerate(Rg):

    a = str( Rg[i] )
    b = str( Ree[i] )
    c = a + " " + b

    file.write(c + '\n')
file.close()
