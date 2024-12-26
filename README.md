# Single Polymer Chain Brownian Dynamics

A Python implementation to simulate single polymer chains using Brownian dynamics. The simulation models beads connected by springs under the influence of Brownian motion, excluded volume effects, and finitely extensible non-linear elastic (FENE) forces.

## Physics Background

The simulation is based on the *overdamped Langevin equation*, which models the motion of particles in a fluid where inertial effects are negligible. The forces acting on each bead include:

- Brownian\random forces (thermal noise)
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
- Configurable number of time steps, (`N = 2000001` by default but can be changed)
- Adjustable simulation parameters (time step, spring constants, etc.)
- Calculates physical properties like:
  - Radius of gyration
  - End-to-end distance
  - Mean square displacement
- Supports different chain configurations:
  - Ideal chains (no excluded volume)
  - Real chains (with excluded volume)
  - Circular chains (polymer rings)


## Usage

### Bead Class

#### Overview
The `Bead` class represents a 2D bead with a radius and an initial position in Cartesian coordinates. It includes methods to initialize the bead and calculate forces acting on it.

#### Attributes
- `x` (float): The x-coordinate of the bead's center.
- `y` (float): The y-coordinate of the bead's center.
- `r` (float): The radius of the bead.

#### Methods

##### `__init__(self, x=0, y=0, r=.04)`
Initializes a new instance of the `Bead` class with the specified position and radius.

**Parameters:**
- `x` (float): Initial x-coordinate of the bead's center. Default is `0`.
- `y` (float): Initial y-coordinate of the bead's center. Default is `0`.
- `r` (float): Radius of the bead. Default is `0.04`.

##### `force_calculate(self, j, jj=None, k=0, k_ev=0, Ls=None, kBT=1, lk=1, conf='linear')`
Computes the forces acting on the bead.

**Parameters:**
- `j` (int): Index of the current bead.
- `jj` (int, optional): Another index. Default is `None`.
- `k` (float): Spring constant. Default is `0`.
- `k_ev` (float): Excluded volume interaction constant. Default is `0`.
- `Ls` (float, optional): Contour length. Default is `None`
- `kBT` (float): Thermal energy. Default is `1`.
- `lk` (float): Kuhn length. Default is `1`.
- `conf` (str): Configuration type, either `'linear'` or `'circular'`. Default is `'linear'`.
  - Please note that this parameter only affects a chain of beads and has no barring on an individual bead.

##### `advance(self, Δt, b=1, κ=0)`
Advances the bead's position over time based on the sum of forces acting on it, following the overdamped Langevin equation.

**Parameters:**
- `Δt` (float): time step size.
- `b` (float): Drag coefficient. Default is `1`.
- `κ` (float): Spring constant (replaces `k` to avoid kernel confusion). Default is `0`.

**Returns:**
None, but updates:
- `positions_xy`: List of (x, y) coordinates over time
- `all_pos_xy`: Global list containing position histories of all beads
- Updates the bead's position attributes (`self.x`, `self.y`)

**Algorithm:**
1. Initializes with current position
2. For each time step:
   - Calculates total force using `force_calculate`
   - Updates position using overdamped Langevin equation
   - Stores new position
3. Converts position history to numpy array

#### Example Usage
```python
# Create a new bead with initial position (0, 0) and radius 0.04
bead = Bead(x=0, y=0, r=0.04)

# Calculate forces on the bead with index 1
bead.force_calculate(j=1)

# Advance the bead's position with time step 0.0001
bead.advance(Δt=0.0001, b=1)
```

### Simulation Class

#### Overview
The `Simulation` class manages the simulation of a Brownian polymer chain consisting of multiple beads. It supports both linear and circular chain configurations.

#### Attributes
- `nbeads` (int): Number of beads in the polymer chain
- `conf` (str): Configuration type (`'linear'` or `'circular'`)
- `ψ` (float): Initial angle between beads (circular configuration only)
- `ρ` (float): Initial radius of circular configuration = `(0.09*nbeads)/(2*π)`
- `beads` (list): List of Bead objects comprising the polymer chain

#### Methods

##### `__init__(self, nbeads, x=0, y=0, conf='linear')`
Initializes a new polymer chain simulation.

**Parameters:**
- `nbeads` (int): Number of beads in the chain
- `x` (float): Initial x-coordinate. Default is `0`.
- `y` (float): Initial y-coordinate. Default is `0`.
- `conf` (str): Configuration type (`'linear'` or `'circular'`). Default is `'linear'`.

**Initialization Process:**
1. Sets up chain configuration
2. Creates global Brownian forces arrays (`Fx_sim`, `Fy_sim`)
3. Initializes beads with appropriate spacing:
   - Linear: Beads spaced `0.09` units apart horizontally
   - Circular: Beads arranged in a circle with radius `ρ`

##### `init_bead(self, x=0, y=0)`
Creates a new bead instance.

**Parameters:**
- `x` (float): x-coordinate for the bead. Default is `0`.
- `y` (float): y-coordinate for the bead. Default is `0`.

**Returns:**
- `Bead`: A new Bead instance

##### `advance(self, Δt, b=1, κ_ev=0)`
Advances the entire polymer chain simulation over time.

**Parameters:**
- `Δt` (float): time step size
- `b` (float): Drag coefficient. Default is `1`.
- `κ_ev` (float): Excluded volume interaction constant. Default is `0`.

**Updates:**
- `xs`, `ys`: Lists of current bead positions
- `end_to_end`: End-to-end distance at each time step
- `Rg`: Radius of gyration at each time step
- `all_sim_pos`: Position history of all beads

**Algorithm:**
1. Initializes position tracking arrays
2. Records initial configuration
3. For each time step:
   - Updates positions of all beads simultaneously
   - Calculates and stores chain properties (end-to-end distance, radius of gyration)
4. Resets position trackers after simulation

#### Example Usage
```python
# Create a linear polymer chain with 10 beads
sim = Simulation(nbeads=10, conf='linear')

# Run simulation with excluded volume interactions
sim.advance(Δt=0.0001, κ_ev=200)

# Create a circular polymer chain
circular_sim = Simulation(nbeads=10, conf='circular')

# Run simulation
circular_sim.advance(Δt=0.0001)
```

#### Global Variables Affected
- `Fx_sim`, `Fy_sim`: Arrays of random forces
- `xs`, `ys`: Current positions of all beads
- `end_to_end`: Chain end-to-end distance history
- `Rg`: Radius of gyration history
- `all_sim_pos`: Complete position history of all beads

#### Notes
- Standard simulation parameters:
  - `Ls = 1` (contour length)
  - `lk = 0.1` (Kuhn length)
  - `kBT = 1` (thermal energy)
  - `N` (number of time steps)
  - `Δt` (time step size)
    - `N` and `Δt` determine the total simulation time and temporal resolution. The smaller `Δt` is:
      - More accurate simulation (better resolution of the dynamics)
      - But requires more steps (larger `N`) to simulate the same total time
- Optimal excluded volume constant (`κ_ev`) is typically around `200`
- Position updates are synchronized across all beads to maintain chain integrity
- The simulation uses the overdamped Langevin equation for bead motion

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

This project is licensed under the MIT License. See the [LICENSE](https://github.com/spyderkam/brownian-motion/blob/main/LICENSE) file for details.
