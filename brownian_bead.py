#!/usr/bin/env python3

# 3/11/2021 — 9/7/2022

import numpy as np

###############################################################################

# preamble

N = 2000001                        # N = number of time steps (iterations)
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

  def __init__(self, x=0, y=0, r=.04):
    self.x = x     # x position of bead's center
    self.y = y     # y position of bead's center
    self.r = r     # bead's radius, best not to play with this

  def force_calculate(self, j, jj=None, k=0, k_ev=0, Ls=None, kBT=1, lk=1, conf='linear'):
    """Compute the forces on each bead."""

    if xs != None:
      if conf == 'circular':
        r1 = np.sqrt( (self.x - xs[j-1])**2 + (self.y - ys[j-1])**2 )
        if j == (len(xs) - 1):
          r2 = np.sqrt( (self.x - xs[0])**2 + (self.y - ys[0])**2 )
          FENE2 = (3*kBT/lk)*( (r2/Ls) / (1 - (r2/Ls)**2) )
          x_FENE2 = FENE2 * (self.x - xs[0])/r2     # FENE2 * cos(θ)
          y_FENE2 = FENE2 * (self.y - ys[0])/r2     # FENE2 * sin(θ)
        else:
          r2 = np.sqrt( (self.x - xs[j+1])**2 + (self.y - ys[j+1])**2 )
          FENE2 = (3*kBT/lk)*( (r2/Ls) / (1 - (r2/Ls)**2) )
          x_FENE2 = FENE2 * (self.x - xs[j+1])/r2     # FENE2 * cos(θ)
          y_FENE2 = FENE2 * (self.y - ys[j+1])/r2     # FENE2 * sin(θ)
        
        FENE1 = (3*kBT/lk)*( (r1/Ls) / (1 - (r1/Ls)**2) )     
        x_FENE1 = FENE1 * (self.x - xs[j-1])/r1     # FENE1 * cos(φ)
        y_FENE1 = FENE1 * (self.y - ys[j-1])/r1     # FENE1 * sin(φ)

        x_force = Fx_sim[jj, j] - (x_FENE1 + x_FENE2)
        y_force = Fy_sim[jj, j] - (y_FENE1 + y_FENE2)

      elif conf == 'linear':
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

      if k_ev:
        xj = np.delete(np.array(xs), j)     # temp list w/all current pos - self.x
        yj = np.delete(np.array(ys), j)
        biscl_x = []     # list of (b)ead (i)nteraction (s)o-(c)alled '(l)engths'
        biscl_y = []

        for index, xpos in enumerate(xj):
          Δx = self.x - xj[index]
          Δy = self.y - yj[index]
          d = np.sqrt( Δx**2 + Δy**2 )

          if d < Ls/1.336:     # check if any volume is excluded
            biscl_x.append((Ls/1.336 - d) * Δx/d)
            biscl_y.append((Ls/1.336 - d) * Δy/d)
              
        x_force = x_force + k_ev*np.sum(biscl_x)
        y_force = y_force + k_ev*np.sum(biscl_y);

    elif xs == None:     # I guess this is just for simple Brownian motion
      x_force = Fx[j] + k*self.x
      y_force = Fy[j] + k*self.y
    return (x_force, y_force)

  def advance(self, Δt, b=1, κ=0):     # κ instead of k just in case the kernel gets confused
    """Advance the beads's position based on ΣF."""

    positions_xy = [(self.x, self.y)]     # initialize
    for i in range(N-1):     # len()-1 cuz already have the initial entry
      self.x = self.x + (self.force_calculate(k=κ, j=i)[0] / b)*Δt          # advance position
      self.y = self.y + (self.force_calculate(k=κ, j=i)[1] / b)*Δt          # x[i] = x[i-1] + (F[i]/b)*Δt

      positions_xy.append( (self.x, self.y) )                               # store the advanced positions
    all_pos_xy.append(positions_xy)           # append current bead pos in all_pos_xy
    positions_xy = np.array(positions_xy)     # positions_xy list is now a NumPy array
##################################################


################  Simulation Class  ################
class Simulation:
  """Basic simulation of Brownian polymer chain. Based on Bead class."""

  def __init__(self, nbeads, x=0, y=0, conf='linear'):
    self.nbeads = nbeads
    self.conf = conf

    if conf == 'circular':
      self.ψ = (2*np.pi)/nbeads     # initial angle between beads
      self.ρ= 0.09/self.ψ           # radius = (0.09*nbeads)/(2*np.pi)
      self.beads = [self.init_bead(self.ρ*np.cos(self.ψ*i), self.ρ*np.sin(self.ψ*i)) for i in range(nbeads)]
    
    elif conf == 'linear':
       self.beads = [self.init_bead(i*.09,0) for i in range(nbeads)]

    global Fx_sim     # global forces for simulation
    global Fy_sim

    Fx_sim = []     # initialize the Brownian forces...
    Fy_sim = []     # ... making it = nbeads x N
    for i in range(N-1):     # N vs N-1
      Fx_sim.append(np.random.normal(0, 100, nbeads))
      Fy_sim.append(np.random.normal(0, 100, nbeads))
    Fx_sim = np.array(Fx_sim)
    Fy_sim = np.array(Fy_sim)

  def init_bead(self, x=0, y=0):
    return Bead(x, y)

  def advance(self, Δt, b=1, κ_ev=0):
    """Advance the simulation."""

    global xs   # not sure why but must globalize to reflect global change
    global ys

    xs = []  # list containing the current x positions of all beads in sim...
    ys = []  # ...all the particles move at once

    xj = []  # store new positions here to avoid changing xs before all beads...
    yj = []  # ... have advanced. Then set xs = xj so all beads advance at once

    for i in range(self.nbeads):
      if self.conf == 'circular':
        all_sim_pos.append([(self.ρ*np.cos(self.ψ*i), self.ρ*np.sin(self.ψ*i))])
      elif self.conf == 'linear':
        all_sim_pos.append([(i*.09,0)]) 

    for bead in self.beads:
      xs.append(bead.x)     # store all the init pos of the beads
      ys.append(bead.y)
    end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))  # 1st e2e element
    Rg.append(np.sqrt(np.var(xs) + np.var(ys)))          # 1st Rg element

    for i in range(N-1):
      for n, bead in enumerate(self.beads):                                                           # Ls=1.5d0, lk=d0
        bead.x = bead.x + (bead.force_calculate(k_ev=κ_ev, j=n, jj=i, Ls=1, lk=.1, kBT=1, conf=self.conf)[0] / b)*Δt  # Ls=1, lk=.1
        bead.y = bead.y + (bead.force_calculate(k_ev=κ_ev, j=n, jj=i, Ls=1, lk=.1, kBT=1, conf=self.conf)[1] / b)*Δt  # Ls=.4 & lk=.04

        xj.append(bead.x); yj.append(bead.y)
        all_sim_pos[n].append( (bead.x, bead.y) )
      xs = xj; ys = yj       # advance all at once
      xj = []; yj = []       # reset to advance all at once next time
      end_to_end.append((xs[-1] - xs[0], ys[-1] - ys[0]))   # e2e @ each time
      Rg.append(np.sqrt(np.var(xs) + np.var(ys)))           # Rg @ each time; Rg^2 would be better for computation time
    xs = None; ys = None     # reset xs and ys to
###############################################################################


if __name__ == '__main__':
  '''The Fundumentals'''

  import matplotlib; matplotlib.use('TkAgg')     # try disabling if error
  import matplotlib.pyplot as plt
  from scipy.optimize import curve_fit

  Δt = .0001      # time step
  t = []          # time array
  tt = 0          # temp element of time array
  for i in range(N):
    t.append(tt)
    tt += Δt

  sim = Simulation(nbeads=2)
  sim.advance(Δt, κ_ev=200)     # κ_ev = 200 seems optimal

  Ree = []     # end-to-end radius; i.e., \sqrt{x^2 + y^2}
  for xy in end_to_end:
    Ree.append(np.sqrt(xy[0]**2 + xy[1]**2))

  '''Some More Fundumentals'''

  def compute_MSD(positions_xy):
    totalsize = len(positions_xy)
    msd = []     # initialize

    for i in range(totalsize):
      j = i + 1

      if totalsize != j:     # don't want a division by zero in class method
        msd.append(np.sum( (positions_xy[j::] - positions_xy[0:-j])**2 ) / (totalsize -j))
    return np.array(msd)

  all_x = []     # list containing all x positions of all beads at each time
  all_y = []
  for i, bead in enumerate(np.array(all_sim_pos)):
    all_x.append([])     # append empty lists equal to the number of beads
    all_y.append([])
    for j, s in enumerate(bead):      # consider all the positions of each individual bead
      all_x[i].append(bead[j, 0])     # append the x pos at each time (i.e., x pos at time j)
      all_y[i].append(bead[j, 1])
  # all_x and all_y lists now complete

  print("no obvious errors...")
