# Brownian Polymer Chain Simulation

A Python implementation to simulate single polymer chains using Brownian dynamics. The simulation models beads connected by springs under the influence of Brownian motion, excluded volume effects, and finitely extensible non-linear elastic (FENE) forces.

## Physics Background

The simulation is based on the overdamped Langevin equation, which models the motion of particles in a fluid where inertial effects are negligible. The forces acting on each bead include:

- Random Brownian forces (thermal noise)
- FENE spring forces between consecutive beads
- Optional excluded volume forces between beads

The key equation governing the bead motion is:


```
dx/dt = ∑F/b
```

where:
- `x` is position
- `F` represents the sum of forces
- `b` is the drag coefficient

## Features

- Simulates both linear and circular polymer chains
- Configurable number of beads
- Adjustable simulation parameters (time step, spring constants, etc.)
- Calculates physical properties like:
  - Radius of gyration
  - End-to-end distance
  - Mean square displacement
- Supports different chain configurations:
  - Ideal chains (no excluded volume)
  - Real chains (with excluded volume)
  - Circular chains (polymer rings)

## Installation

1. Clone this repository:
```bash
git clone [repository-url]
```

2. Required dependencies:
```bash
pip install numpy matplotlib
```

## Usage

The main classes are:

### Bead Class
```python
bead = Bead(x=0, y=0, r=0.04)  # Create a bead at (0,0) with radius 0.04
```

### Simulation Class
```python
sim = Simulation(nbeads=10, conf='linear')  # Create a linear chain with 10 beads
sim.advance(dt=0.0001, κ_ev=200)  # Run simulation with excluded volume
```

### Example
```python
import brownian_bead as bb

# Create and run a simple simulation
sim = bb.Simulation(nbeads=2)
sim.advance(dt=0.0001, κ_ev=200)

# Access simulation results
print(bb.Rg)  # Radius of gyration
print(bb.end_to_end)  # End-to-end distances
```

## Parameters

Key simulation parameters:

- `dt`: Time step size (default: 0.0001)
- `N`: Number of iterations (default: 2000001)
- `κ_ev`: Excluded volume constant (default: 200)
- `b`: Drag coefficient (default: 1)
- `Ls`: Maximum spring length (default: 1)
- `lk`: Kuhn length (default: 0.1)
- `kBT`: Temperature energy scale (default: 1)

## Output

The simulation outputs several physical quantities:

- Time series of bead positions
- Radius of gyration (`Rg`)
- End-to-end distance
- Chain conformations at different time points

## Limitations

- 2D simulations only
- No hydrodynamic interactions
- Limited to single chain dynamics

## Contributing

Feel free to open issues or submit pull requests for:

- New features
- Documentation improvements
- Performance optimizations

## License

This project is licensed under the MIT License. See the (LICENSE)[https://github.com/spyderkam/brownian-motion/blob/main/LICENSE] file for details.
