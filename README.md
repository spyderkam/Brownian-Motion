# README (readme needs updating)

Bead class to define Brownian particles. Simulation class based on Bead class for multi-Bead simulations.

<!--- * `brownian_bead_FENE` is also a simulation for the ideal chain/Rouse model but for a FENE force spring (see Wiki); which is what is ultimately desired. --->
* `brownian_bead.py` is the code which simulates a polymer chain with Brownian forces, FENE force, and optional excluded volume force.

  * `run.py` executes `brownian_bead` simulation(s) and saves the coordinates, radius of gyration, and end-to-end radius of the polymer chain.

  * `hist.py` plots histogram of the angle distribution of the chain at different points in time.
