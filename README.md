# README (readme needs updating)

Bead class to define Brownian particles. Simulation class based on Bead class for multi-Bead simulations.

<!--- * `brownian_bead_FENE` is also a simulation for the ideal chain/Rouse model but for a FENE force spring (see Wiki); which is what is ultimately desired. --->
  * `brownian_bead.py` is the code which simulates a linear polymer chain with Brownian forces, FENE force, and optional excluded volume force; no HI force.

  * `run.py` executes `brownian_bead` simulation(s) and saves the coordinates, radius of gyration, and end-to-end radius of the polymer chain.

  * `hist.py` plots histogram of the angle distribution of the chain at different points in time.
  
  
  #### TCPP
  `TCPP.pdf` is my (unofficial) master's thesis which is based on my professor's program for a discrete worm-like chain, the `rpbd` repository, and largely this repository.
