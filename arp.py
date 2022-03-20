# Due to the 'awk' command passed to subprocess.call method, this script can only be executed using Linux/Unix operating systems.

from rpbd import *
import subprocess     # https://www.kite.com/python/answers/how-to-run-bash-commands-in-python
# https://stackoverflow.com/questions/30928540/append-to-file-using-subprocess-in-python
   
m = 14    # number of beads per ring
n = 2      # number of rings
setwrite = 1; run = 100000; N = int(run/setwrite)     # copy/paste

print("now running simulation...")
positions = np.zeros((n*m,3))
radius = 4.0*float(m)/(2.0*np.pi)
for i in range(n):
    for j in range(m):
        positions[i*m+j,0] = 2.0*i
        positions[i*m+j,1] = radius*np.cos(2.0*np.pi*j/float(4))
        positions[i*m+j,2] = radius*np.sin(2.0*np.pi*j/float(4))
    
# create simulation object and set positions
sim = Simulation(numRings=n, beadsPerRing=m, hiType=0)
sim.setPositions(positions)
sim.setOutput("simulated_data/No-HI/inst-temp.xyz")
sim.setCalc(1)
sim.setWrite(setwrite)
result = sim.run(run)

sim.prepOutput()
sim.integrateInst(run)     # sim.integratePre(100000) 
sim.closeOutput()


# running bash command
subprocess.call(f"awk '!x[$0]++' simulated_data/No-HI/inst-temp.xyz > simulated_data/No-HI/temp-output_n={n}_m={m}.txt", shell=True)


coos = np.genfromtxt(f'simulated_data/No-HI/temp-output_n={n}_m={m}.txt', skip_header=2)  # list of coordinates
#print(len(coos), len(coos)//(m*n), m*n, 3)  #
coos = coos.reshape(N, m*n, 3)     # len(coos)//(m*n)       

# radius of gyration
Rg = []                   # initialize
var_list = [[],[],[]]     # variance list
xyz_list = [[],[],[]]     # list of x,y,z
for time_step in coos:
    for bead in time_step:
        for i, coo in enumerate(bead):
            xyz_list[i].append(coo)
    
    for i, xyz in enumerate(xyz_list):
        var_list[i].append(np.var(xyz))
        
    Rg.append(np.sqrt(np.sum(var_list)))     # append Rg for current time step iteration
              
    var_list = [[],[],[]]     # reset
    xyz_list = [[],[],[]]     # reset   
    
print("...now writing to file")   
file = open(f"simulated_data/No-HI/ave_rg_nm{n}{m}.dat", "a")

βmax = [np.array(Rg, dtype=object), np.array(coos[:,:,0], dtype=object), np.array(coos[:,:,1], dtype=object), 
        np.array(coos[:,:,2], dtype=object)]     # arange data for output, didn't know what else to call it...

for i, s in enumerate(Rg):  
    
    a = str( βmax[0][:][i] ).replace("\n", "")
    b = str( βmax[1][:][i] ).replace("\n", "")
    c = str( βmax[2][:][i] ).replace("\n", "")
    d = str( βmax[3][:][i] ).replace("\n", "")
    
    e = a + "\t" + b + "\t" + c + "\t" + d
    e = e.replace("[","").replace("]","")
    
    file.writelines(e + "\n")
file.close()
