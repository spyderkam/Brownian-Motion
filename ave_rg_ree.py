from brownian_bead import *

Δt = 0.001
t = []
tt = 0
for i in range(N):
    t.append(tt)
    tt += Δt

print("now running simulation...")
    
sim0 = Simulation(nbeads=100)
sim0.advance(Δt, κ_ev=1) 

Ree = []
for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

print("...now writing to file")   

file = open("simulated_data/ave_rg_ree_n100.dat", "w")
for i, s in enumerate(Rg):

    a = str( Rg[i] )
    b = str( Ree[i] )
    c = a + " " + b

    file.write(c + '\n')
file.close()