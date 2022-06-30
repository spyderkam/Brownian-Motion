#!/usr/bim/env python3

from brownian_bead import *
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#import matplotlib as mpl
#matplotlib.rcParams['lines.markersize'] = 150     # 140 good for nbeads=2

Δt = .001
beads = 7
sim = Simulation(nbeads=beads)
sim.advance(Δt, κ_ev=3)

arr = []
for i in range(N):
    arr.append([])
for i, bead in enumerate(all_sim_pos):
    for j, time_point in enumerate(bead):
        arr[j].append(all_sim_pos[i][j])
positions = np.array(arr)

def init():
    scatterplot.set_offsets([[], []])
    return [scatterplot]

def update(i, scatterplot, positions):
    scatterplot.set_offsets(positions[i])
    return [scatterplot]

fig = plt.figure()

plt.xlim(-.01, (beads-1)*.09 + .01)     # -.06, .14 good for 2 beads
plt.ylim(-.2, .2)                     # -.75, .75 good for 2 beads

#points_whole_ax = 5 * .8 * 72    # 1 point = dpi / 72 pixels
#radius = .04
#points_radius = 2 * radius / 1 * points_whole_ax     # 1*"the ratio of decrease from 1"

scatterplot = plt.scatter([], [], marker='o', color='blue')

anim = animation.FuncAnimation(fig, update, init_func=init, fargs=(scatterplot, positions), interval=200, frames=int(N/10), blit=True, repeat=True)
#plt.grid()
plt.axis('off')
f = r"animationpy.mov" 
writervideo = animation.FFMpegWriter(fps=60) 
anim.save(f, writer=writervideo)
plt.show()



"""
Recources
[1] https://stackoverflow.com/questions/49835134/matplotlib-animate-multiple-scatter-plots
[2] https://stackoverflow.com/questions/33094509/correct-sizing-of-markers-in-scatter-plot-to-a-radius-r-in-matplotlib
"""
