# Brownian Motion Polymer Chain Simulation

A Python-based simulation of Brownian motion for polymer chains, focusing on studying chain dynamics through metrics like radius of gyration (Rg) and end-to-end distance (Ree).

## Overview

This project simulates the behavior of polymer chains using Brownian dynamics, allowing for both linear and circular chain configurations. The simulation accounts for various physical interactions including FENE (Finitely Extensible Nonlinear Elastic) potential and volume exclusion.

## Core Features

- Multiple chain configurations supported:
  - Linear chains
  - Circular chains
- Configurable parameters:
  - Number of beads (tested from $5$ to $100$ beads)
  - Time steps (default $N = 2,\\\!000,\\\!001$)
  - Time step size (Δt = 0.0001)
  - Volume exclusion strength (`k_ev`, optimal at `200`)
- Analysis capabilities:
  - Radius of gyration ($R_{\mathrm{g}}$) calculation
  - End-to-end distance ($R_{\mathrm{ee}}$) measurement
  - Mean Square Displacement ($\mathrm{MSD}$) computation
  - Relaxation time ($\tau$) extraction

## Project Structure

```
project/
│
├── brownian_bead.py    # Core simulation engine
├── run.py              # Script for running simulations and saving data
└── auto-extract.ipynb  # Analysis notebook for processing results
```

### File Descriptions

#### brownian_bead.py
- Main simulation engine
- Contains `Bead` and `Simulation` classes
- Implements force calculations and position updates
- Handles both linear and circular chain configurations

#### run.py
- Executes simulations for specified chain lengths
- Saves position data and metrics ($R_{\mathrm{g}}$ & $R_{\mathrm{ee}}$) to files
- Supports batch processing of multiple chain lengths

#### auto-extract.ipynb
- Jupyter notebook for data analysis
- Extracts relaxation times ($\tau$)
- Calculates average $R_{\mathrm{g}}$ and $R_{\mathrm{ee}}$ after relaxation
- Includes various data processing options (e.g., 85/15 split for analysis)

## Usage

1. Set up simulation parameters in `brownian_bead.py`:
```python
N = 2000001  # Number of time steps
Δt = 0.0001  # Time step size
```

2. Run simulation using `run.py`:
```python
python run.py
```

3. Analyze results using `auto-extract.ipynb`

## Implementation Details

### Force Calculations
- FENE potential for chain connectivity
- Volume exclusion interactions between beads
- Random Brownian forces in x and y directions

### Data Storage
- Position data saved in numpy format
- Metrics (Rg & Ree) saved in text files
- Organized directory structure for different configurations

### Analysis Methods
- Exponential decay fitting for relaxation times
- Statistical averaging of steady-state properties
- Various time window options for data analysis (5/95, 15/85, 25/75 splits)

## Requirements

- Python 3.10 or higher
- NumPy
- Matplotlib
- SciPy
- Jupyter (for analysis notebook)
- TkAgg backend for matplotlib (configurable)

## Notes

- Default volume exclusion strength (k_ev = 200) has been optimized through testing
- The simulation supports both single and multiple trial runs
- Data analysis typically focuses on the final 85% of simulation time
- Position data can be quite large due to the number of time steps

## Author

spyderkam (March 11, 2021 - September 7, 2022)
