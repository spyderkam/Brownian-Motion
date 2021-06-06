#!/usr/bin/env python3

import numpy as np
import random

###############################################################################

# preamble

parity = (-1, +1)
N = 100001                               # N = number of steps
Fx = np.random.normal(0.5, 0.2, N)     # random forces in the x-direction...
Fy = np.random.normal(0.5, 0.2, N)     # ...for the Bead class   ...   †

all_pos_xy = []     # (x, y) of bead in random walk
all_vel_xy = []     # (vx, vy) of bead in random walk
all_sim_pos = []    # (x, y) position of all beads in simulation
#all_sim_vel = []
xs = None; ys = None
end_to_end = []     # end 2 end dist. of chain @ each time in simulation
Rg = []             # radius of gyration


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


  def force_calculate(self, j, jj=None, k=0):     # avoid calling outside of 'advance' method
    if xs != None:  #new loop
      if j == 0:
        x_force = Fx_sim[jj, j] - k*(self.x - xs[j+1])  # 1 vs j+1
        y_force = Fy_sim[jj, j] - k*(self.y - ys[j+1])
      elif j == (len(xs) - 1):
        x_force = Fx_sim[jj, j] - k*(self.x - xs[j-1])
        y_force = Fy_sim[jj, j] - k*(self.y - ys[j-1])
      else: #if j != (0 or len(xs)-1):
        x_force = Fx_sim[jj, j] * random.choice(parity) - k*(2*self.x - xs[j-1] - xs[j+1])
        y_force = Fy_sim[jj, j] * random.choice(parity) - k*(2*self.y - ys[j-1] - ys[j+1])
      #x_force = Fx[j] * random.choice(parity) - k*(len(xs)*self.x - np.sum(xs))  #old #j=i
      #y_force = Fy[j] * random.choice(parity) - k*(len(ys)*self.y - np.sum(ys))  #old

    elif xs == None:
      x_force = Fx[j] * random.choice(parity) - k*self.x     # parity is random EACH time
      y_force = Fy[j] * random.choice(parity) - k*self.y     # j - 1 vs j ???   ...   †parity no need...
    return (x_force, y_force)


  def advance(self, Δt, b=1, κ=0):     # κ instead of k just in case the kernel gets confused
    """Advance the beads's position according to its velocity...
    ...need function to change vs based on F."""

    positions_xy = [(self.x, self.y)]     # initialize
    velocities = [(self.vx, self.vy)]     # initialize

    for i in range(N-1):     # len()-1 cuz already have the initial entry
      self.x = self.x + (self.force_calculate(k=κ, j=i)[0] / b)*Δt          # advance position
      self.y = self.y + (self.force_calculate(k=κ, j=i)[1] / b)*Δt          # x[i] = x[i-1] + (F[i]/b)*Δt

      positions_xy.append( (self.x, self.y) )                               # store the advanced positions
      velocities.append( (self.x * Δt + self.vx, self.y*Δt + self.vy) )     # store the advanced velocities

      self.vx = np.array(velocities)[i, 0]                                  # advance velocity based...
      self.vy = np.array(velocities)[i, 1]                                  # ...: vy[i] = y[i]*Δt + vy[i-1]

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
    self.beads = [self.init_bead(x=i/10, y=0) for i in range(nbeads)]     # i vs 0

    global Fx_sim     # global forces for simulation
    global Fy_sim

    Fx_sim = []
    Fy_sim = []
    for i in range(N):
      Fx_sim.append(np.random.normal(0, 1, nbeads))     # SCALE W/κ
      Fy_sim.append(np.random.normal(0, 1, nbeads))     # SCALE W/κ
    Fx_sim = np.array(Fx_sim)
    Fy_sim = np.array(Fy_sim)

  def init_bead(self, x=0, y=0, vx=0, vy=0):
    return Bead(x, y, vx, vy)


  def advance(self, Δt, b=1, κ=0):
    global xs   # not sure why but must globalize to reflect global change
    global ys

    xs = []  # list containing the current x positions of all beads in sim...
    ys = []  # ...all the particles move at once

    xj = []  # Store new positions here to avoid changing xs before all beads...
    yj = []  # ... have advanced. Then set xs = xj so all beads advance at once

    for i in range(self.nbeads):
      all_sim_pos.append([(i/10,0)])     # [(0,0)] vs [] vs [(i,i)]
      #all_sim_vel.append([])

    for bead in self.beads:
      xs.append(bead.x)     # store all the init pos of the beads
      ys.append(bead.y)
    end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))  # 1st e2e element
    Rg.append(np.sqrt(np.var(xs) + np.var(ys)))          # 1st Rg element

    for i in range(N-1):
      for n, bead in enumerate(self.beads):
        bead.x = bead.x + (bead.force_calculate(k=κ, j=n, jj=i)[0] / b)*Δt
        bead.y = bead.y + (bead.force_calculate(k=κ, j=n, jj=i)[1] / b)*Δt
        #xs.append(bead.x); ys.append(bead.y)     # not all at once...
        #xs = xs[1:]; ys = ys[1:]                 # ...delete one by one
        xj.append(bead.x); yj.append(bead.y)
        all_sim_pos[n].append( (bead.x, bead.y) )
        #all_sim_vel[n].append( () )     #  might eventually include velocity
      xs = xj; ys = yj       # advance all at once
      xj = []; yj = []       # reset to advance all at once next time
      end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))   # e2e @ each time
      Rg.append(np.sqrt(np.var(xs) + np.var(ys)))           # Rg @ each time
    xs = None; ys = None     # reset xs and ys to
###############################################################################


if __name__ == '__main__':
  # The Fundumentals
  import matplotlib; matplotlib.use('TkAgg')
  import matplotlib.pyplot as plt
  from scipy.optimize import curve_fit

  Δt = 0.01
  t = []
  tt = 0
  for i in range(N):
    t.append(tt)
    tt += Δt

  sim0 = Simulation(nbeads=2)
  sim0.advance(Δt, κ=1)

  Ree = []     # end to end radius; i.e., \sqrt{x^2 + y^2}
  for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))
