from brownian_bead import *

nbeads = 25

Δt = .0001
t = []
tt = 0
for i in range(N):
    t.append(tt)
    tt += Δt

print(f"now running simulation for {nbeads} beads, {N} time-steps, each of size {Δt}...")

sim = Simulation(nbeads=nbeads)
sim.advance(Δt, κ_ev=0.0)


###########################
print("...now saving positions...")

array0 = np.array(all_sim_pos)
with open(f'path/to/nbeads{NoB}.npy', 'wb') as f:
    np.save(f, array0)
array1 = np.load(f'path/to/nbeads{NoB}.npy')

# make sure arrays are the same
comparison = (array0 == array1)     # don't need parantheses
equal_arrays = comparison.all()
if equal_arrays != True:
    print(equal_arrays)
###########################


Ree = []
for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

print("...now writing Rg and Ree to file")

file = open(f"path/to/ave_rg_ree_n{NoB}.dat", "w")
for i, s in enumerate(Rg):

    a = str( Rg[i] )
    b = str( Ree[i] )
    c = a + " " + b

    file.write(c + '\n')
file.close()
