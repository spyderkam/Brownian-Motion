#!/usr/bin/env python3

# 6/5/2021 — 7/7/2021

import numpy as np

###############################################################################

# preamble

N = 100001                         # N = number of time steps
Fx = np.random.normal(0, 1, N)     # random forces in the x-direction...
Fy = np.random.normal(0, 1, N)     # ...for the Bead class

all_pos_xy = []     # (x, y) of bead in random walk
all_sim_pos = []    # (x, y) position of all beads in simulation
xs = None; ys = None
end_to_end = []     # end 2 end dist. of chain @ each time in simulation
Rg = []             # radius of gyration


################  Bead Class  ################
class Bead:
  """A 2D bead w/radius and initial position in Cartesian coordinates."""

  def __init__(self, x=0, y=0, r=.08):
    self.x = x     # x position of bead's center
    self.y = y     # y position of bead's center
    self.pos = (self.x, self.y)
    self.r = r     # bead radius



  def force_calculate(self, j, jj=None, k=0, k_ev=0, Ls=None, kBT=1, lk=1):
    if xs != None:
      if j == 0:
        r = np.sqrt( (self.x - xs[j+1])**2 + (self.y - ys[j+1])**2 )
        FENE = (3*kBT/lk)*( (r/Ls) / (1 - (r/Ls)**2) )
        x_FENE = FENE * (self.x - xs[j+1])/r     # FENE * cos(θ)
        y_FENE = FENE * (self.y - ys[j+1])/r     # FENE * sin(θ)

        x_force = Fx_sim[jj, j] - x_FENE
        y_force = Fy_sim[jj, j] - y_FENE
      elif j == (len(xs) - 1):
        r = np.sqrt( (self.x - xs[j-1])**2 + (self.y - ys[j-1])**2 )
        FENE = (3*kBT/lk)*( (r/Ls) / (1 - (r/Ls)**2) )
        x_FENE = FENE * (self.x - xs[j-1])/r     # FENE * cos(φ)
        y_FENE = FENE * (self.y - ys[j-1])/r     # FENE * sin(φ)

        x_force = Fx_sim[jj, j] - x_FENE
        y_force = Fy_sim[jj, j] - y_FENE
      else:
        r1 = np.sqrt( (self.x - xs[j-1])**2 + (self.y - ys[j-1])**2 )
        r2 = np.sqrt( (self.x - xs[j+1])**2 + (self.y - ys[j+1])**2 )

        FENE1 = (3*kBT/lk)*( (r1/Ls) / (1 - (r1/Ls)**2) )
        FENE2 = (3*kBT/lk)*( (r2/Ls) / (1 - (r2/Ls)**2) )

        x_FENE1 = FENE1 * (self.x - xs[j-1])/r1     # FENE1 * cos(φ)
        y_FENE1 = FENE1 * (self.y - ys[j-1])/r1     # FENE1 * sin(φ)
        x_FENE2 = FENE2 * (self.x - xs[j+1])/r2     # FENE2 * cos(θ)
        y_FENE2 = FENE2 * (self.y - ys[j+1])/r2     # FENE2 * sin(θ)

        x_force = Fx_sim[jj, j] - (x_FENE1 + x_FENE2)
        y_force = Fy_sim[jj, j] - (y_FENE1 + y_FENE2)

      # new stuff
      xj = list(set(xs) - set([self.x]))     # here xj is the temp list which...
      yj = list(set(ys) - set([self.y]))     #... has all pos apart from self.x
      biscl_x = []     # list of (b)ead (i)nteraction (s)o-(c)alled '(l)engths'
      biscl_y = []
      for index, xpos in enumerate(xj):
        if abs(self.x - xj[index]) < 2*self.r:     # check if volume is excluded
          biscl_x.append(abs(self.x - xj[index]))  # append if it is
        #if abs(self.y - yj[index]) < 2*self.r:    ## not sure why these two...
        #  biscl_y.append(abs(self.y - yj[index])) ## ...lines give error  :/
      for ypos in yj:
        if abs(self.y - ypos) < 2*self.r:
          biscl_y.append(abs(self.y - ypos))

      x_force = x_force + k_ev*np.sum(biscl_x)     # - vs +
      y_force = y_force + k_ev*np.sum(biscl_y)
      ###########

    elif xs == None:     # I guess this is just for simple Brownian motion
      x_force = Fx[j] - k*self.x
      y_force = Fy[j] - k*self.y
    return (x_force, y_force)


  def advance(self, Δt, b=1, κ=0, κ_ev=0):     # κ instead of k just in case the kernel gets confused
    """Advance the beads's position based on ΣF."""

    positions_xy = [(self.x, self.y)]     # initialize

    for i in range(N-1):     # len()-1 cuz already have the initial entry
      self.x = self.x + (self.force_calculate(k=κ, k_ev=κ_ev, j=i)[0] / b)*Δt          # advance position
      self.y = self.y + (self.force_calculate(k=κ, k_ev=κ_ev, j=i)[1] / b)*Δt          # x[i] = x[i-1] + (F[i]/b)*Δt

      positions_xy.append( (self.x, self.y) )                               # store the advanced positions
    all_pos_xy.append(positions_xy)           # append current bead pos in all_pos_xy
    positions_xy = np.array(positions_xy)     # positions_xy list is now a NumPy array
##################################################


################  Simulation Class  ################
class Simulation:
  """Simulation class based on Bead class."""
  def __init__(self, nbeads, x=0, y=0, vx=0, vy=0):
    self.nbeads = nbeads
    self.beads = [self.init_bead(x=i*.09, y=0) for i in range(nbeads)]     # i vs 0

    global Fx_sim     # global forces for simulation
    global Fy_sim

    Fx_sim = []
    Fy_sim = []
    for i in range(N-1):     # N vs N-1
      Fx_sim.append(np.random.normal(0, 1, nbeads))
      Fy_sim.append(np.random.normal(0, 1, nbeads))
    Fx_sim = np.array(Fx_sim)
    Fy_sim = np.array(Fy_sim)

  def init_bead(self, x=0, y=0):
    return Bead(x, y)


  def advance(self, Δt, b=1, κ=0, κ_ev=0):
    """Advance the simulation."""
    global xs   # not sure why but must globalize to reflect global change
    global ys

    xs = []  # list containing the current x positions of all beads in sim...
    ys = []  # ...all the particles move at once

    xj = []  # store new positions here to avoid changing xs before all beads...
    yj = []  # ... have advanced. Then set xs = xj so all beads advance at once

    for i in range(self.nbeads):
      all_sim_pos.append([(i*.09,0)])     # [(0,0)] vs [] vs [(i,i)]

    for bead in self.beads:
      xs.append(bead.x)     # store all the init pos of the beads
      ys.append(bead.y)
    end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))  # 1st e2e element
    Rg.append(np.sqrt(np.var(xs) + np.var(ys)))          # 1st Rg element

    for i in range(N-1):
      for n, bead in enumerate(self.beads):
        bead.x = bead.x + (bead.force_calculate(k=κ, k_ev=κ_ev, j=n, jj=i, Ls=.5, lk=.1)[0] / b)*Δt
        bead.y = bead.y + (bead.force_calculate(k=κ, k_ev=κ_ev, j=n, jj=i, Ls=.5, lk=.1)[1] / b)*Δt

        xj.append(bead.x); yj.append(bead.y)
        all_sim_pos[n].append( (bead.x, bead.y) )
      xs = xj; ys = yj       # advance all at once
      xj = []; yj = []       # reset to advance all at once next time
      end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))   # e2e @ each time
      Rg.append(np.sqrt(np.var(xs) + np.var(ys)))           # Rg @ each time
    xs = None; ys = None     # reset xs and ys to
###############################################################################


if __name__ == '__main__':
  '''The Fundumentals'''
  import matplotlib; matplotlib.use('TkAgg')
  import matplotlib.pyplot as plt
  from scipy.optimize import curve_fit


  def compute_MSD(positions_xy):
    totalsize = len(positions_xy)
    msd = []     # initialize

    for i in range(totalsize):
      j = i + 1

      if totalsize != j:     # don't want a division by zero in class method
        msd.append(np.sum( (positions_xy[j::] - positions_xy[0:-j])**2 ) / (totalsize -j))
    return np.array(msd)


  Δt = 0.01     # time step
  t = []        # time array
  tt = 0        # element of time array
  for i in range(N):
    t.append(tt)
    tt += Δt

  sim0 = Simulation(nbeads=2)
  sim0.advance(Δt, κ=1, κ_ev=5)     # make κ_ev very big and play with it

  Ree = []     # end to end radius; i.e., \sqrt{x^2 + y^2}
  for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

  print("no obvious errors...")
