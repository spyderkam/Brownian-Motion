
# Brownian Motion Polymer Chain Simulation

A Python-based simulation of polymer chains using Brownian dynamics. This project simulates both linear and circular polymer chains with excluded volume interactions.

## Features

- 2D Brownian motion simulation of bead-spring polymer chains
- Supports both linear and circular chain configurations
- Calculates radius of gyration (Rg) and end-to-end distance (Ree)
- Includes excluded volume interactions
- Visualization tools for analyzing chain dynamics

## Files

- `brownian_bead.py`: Core simulation engine with Bead and Simulation classes
- `run.py`: Script to execute simulations and save position/property data
- `hist.py`: Generates histograms of bond angle distributions
- `auto-extract.ipynb`: Jupyter notebook for analyzing Rg and Ree relaxation times

## Usage

1. Run a basic simulation:
```python
import brownian_bead as bb

sim = bb.Simulation(nbeads=20, conf='linear')  # or 'circular'
sim.advance(Δt=0.0001, κ_ev=200)
```

2. Execute predefined simulations:
```bash
python3 run.py
```

3. Generate angle distribution histograms:
```bash
python3 hist.py
```

## Parameters

- `nbeads`: Number of beads in the polymer chain
- `Δt`: Time step for integration
- `κ_ev`: Excluded volume interaction strength
- `conf`: Chain configuration ('linear' or 'circular')

## Data Output

- Position data: Saved as .npy files
- Rg & Ree data: Saved as .dat files
- Histograms: Saved as .png files

## Dependencies

- NumPy
- Matplotlib
- Jupyter (for analysis notebooks)
