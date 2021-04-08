#!/usr/bin/env python3

import numpy as np
import random

###############################################################################


# preamble

parity = (-1, +1)
N = 1001                               # N = number of steps
Fx = np.random.normal(0.5, 0.2, N)     # random forces in the x-direction
Fy = np.random.normal(0.5, 0.2, N)
all_pos_xy = []     # (x, y) of all beads in simulation
all_vel_xy = []     # (vx, vy) of all beads
all_sim_pos = []    # (x, y) position of all beads in simulation
#all_sim_vel = []
xs = None; ys = None #old


################  Bead Class  ################
class Bead:
  """A 2D bead w/initial position and velocity in Cartesian coordinates."""

  def __init__(self, x=0, y=0, vx=0, vy=0):
    self.x = x
    self.y = y
    self.pos = (self.x, self.y)
    self.vx = vx
    self.vy = vy
    self.vel = (vx, vy)



  def force_calculate(self, j, k=0):     # avoid calling outside of 'advance' method
    if xs != None:
      x_force = Fx[j] * random.choice(parity) - k*(len(xs)*self.x - np.sum(xs))
      y_force = Fy[j] * random.choice(parity) - k*(len(ys)*self.y - np.sum(ys))

    elif xs == None:
      x_force = Fx[j] * random.choice(parity) - k*self.x        # parity is random EACH time
      y_force = Fy[j] * random.choice(parity) - k*self.y
    return (x_force, y_force)


  def advance(self, Δt, b=1, κ=0):     # κ instead of k just in case the kernel gets confused
    """Advance the beads's position according to its velocity...
    ...need function to change vs based on F."""

    global velocities         # global to evaluate globally
    global positions          # at this point don't need global for any of these...
    global positions_xy       # ... just nice to have

    positions_xy = [(self.x, self.y)]     # initialize
    velocities = [(self.vx, self.vy)]     # initialize

    for i in range(N-1):     # len()-1 cuz already have the initial entry
      self.x = self.x + (self.force_calculate(k=κ, j=i)[0] / b)*Δt          # advance position
      self.y = self.y + (self.force_calculate(k=κ, j=i)[1] / b)*Δt          # x[i] = x[i-1] + (F[i]/b)*Δt

      positions_xy.append( (self.x, self.y) )                               # store the advanced positions
      velocities.append( (self.x * Δt + self.vx, self.y*Δt + self.vy) )     # store the advanced velocities

      self.vx = np.array(velocities)[i, 0]                                  # advance velocity based...
      self.vy = np.array(velocities)[i, 1]                                  # ...: vy[i] = y[i]*Δt + vy[i-1]

    #positions = []     # initialize in matrix format
    #for position in positions_xy:
    #  r = np.sqrt(position[0]**2 + position[1]**2)  # r = \sqrt{x^2 + y^2}
    #  positions.append(r)                           # list of r at each point
    #positions = np.array(positions)                 # positions list is now a NumPy array

    all_pos_xy.append(positions_xy)     # append current bead pos in all_pos_xy
    all_vel_xy.append(velocities)       # append current bead vel in all_vel_xy

    positions_xy = np.array(positions_xy)     # positions_xy list is now a NumPy array
    velocities = np.array(velocities)         # velocites list is now a NumPy array


  def compute_MSD(self, positions_xy):     # advance first... BEST to def outside class
    totalsize = len(positions_xy)
    msd = []     # initialize

    for i in range(totalsize):
      j = i + 1

      if totalsize != j:     # don't want a division by zero in class method
        msd.append(np.sum( (positions_xy[j::] - positions_xy[0:-j])**2 ) / (totalsize -j))
    return np.array(msd)
##################################################


################  Simulation Class  ################
class Simulation:
  """Simulation class based on Bead class."""
  def __init__(self, nbeads, x=0, y=0, vx=0, vy=0):
    self.nbeads = nbeads
    self.beads = [self.init_bead() for i in range(nbeads)]


  def init_bead(self, x=0, y=0, vx=0, vy=0):
    return Bead(x, y, vx, vy)


  def advance(self, Δt, b=1, κ=0):
    global xs #new  # not sure why but must globalize to reflect global change
    global ys #new

    xs = []  # list containing the current x positions of all beads in sim...
    ys = []  # ...all the particles move at once

    xj = []  # Store new positions here to avoid changing xs before all beads...
    yj = []  # ... have advanced. Then set xs = xj so all beads advance at once

    for i in range(self.nbeads):
      all_sim_pos.append([(0,0)])  # [(0,0)] vs []
      #all_sim_vel.append([])

    for bead in self.beads:
      xs.append(bead.x)     # store all the init pos of the beads
      ys.append(bead.y)

    for i in range(N-1):
      for n, bead in enumerate(self.beads):
        bead.x = bead.x + (bead.force_calculate(k=κ, j=i)[0] / b)*Δt
        bead.y = bead.y + (bead.force_calculate(k=κ, j=i)[1] / b)*Δt
        #xs.append(bead.x); ys.append(bead.y)     # not all at once
        #xs = xs[1:]; ys = ys[1:]                 # delete one by one
        xj.append(bead.x); yj.append(bead.y)
        all_sim_pos[n].append( (bead.x, bead.y) )
        #all_sim_vel[n].append( () )     #  might eventually include velocity
      xs = xj; ys = yj       # advance all at once
      xj = []; yj = []       # reset to advance all at once next time
    xs = None; ys = None     # reset xs and ys to
###############################################################################


if __name__ == '__main__':
  import matplotlib; matplotlib.use('TkAgg')
  import matplotlib.pyplot as plt
  #from scipy.optimize import curve_fit

  time = np.linspace(0, 100, N)
  Δt = time[1] - time[0]

  sim0 = Simulation(nbeads=2)
  sim0.advance(Δt, κ=1)
  b1 = np.array(all_sim_pos[0])
  b2 = np.array(all_sim_pos[1])

  plt.plot(b1[:, 0], b1[:, 1], ',')
  plt.plot(b2[:, 0], b2[:, 1], ',')
  plt.show()