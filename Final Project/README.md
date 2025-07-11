# Final Project: Advanced Vector Fitting with Passivity Enforcement

## Overview
This project combines fast-relaxed vector fitting algorithms with passivity enforcement to create stable, accurate circuit models from S-parameter data. It represents the culmination of the vector fitting and circuit synthesis techniques developed in the previous projects.

## Description
The final project implements a complete workflow from S-parameter data to simulation-ready SPICE netlists. It uses the fast-relaxed vector fitting algorithm to create high-order models, then enforces passivity to ensure time-domain stability. The resulting models are suitable for simulation in Advanced Design Systems (ADS).

## Key Features
- **Fast-Relaxed Vector Fitting**: Advanced fitting algorithm for improved accuracy
- **Passivity Enforcement**: Ensures time-domain stability
- **Multi-port Support**: Handles arbitrary N-port networks
- **SPICE Netlist Generation**: Creates simulation-ready circuits
- **Real-time Visualization**: Shows fitting progress and results
- **Comprehensive Testing**: Includes unit tests and validation

## Core Components
- `Final_Project.py` - Main orchestration script
- `VFDriver.py` - Vector fitting driver and utilities
- `RPDriver.py` - Pole-residue processing driver
- `vectfit3.py` - Core vector fitting algorithm
- `FRPY.py` - Fast-relaxed vector fitting implementation
- `intercheig.py` - Eigenvalue analysis utilities
- `fitcalcABCDE.py`, `fitcalcPRE.py` - Fitting calculation modules
- `rot.py` - Rotation and transformation utilities
- `pr2ss.py` - Pole-residue to state-space conversion
- `violextremaY.py` - Passivity violation detection

## Files
- `Final_Project.py` - Main execution script
- `ex2_Y.mat` - Input S-parameter data
- `data.npz`, `data.np.npz` - Processed data files
- `spice_net_list.sp` - Generated SPICE netlist
- `test/` - Unit test directory
- Various driver and utility modules

## Usage
1. Ensure all required files are in the working directory
2. Run `Final_Project.py`
3. The program will display real-time fitting progress
4. Results include fitted data, error metrics, and SPICE netlist
5. Import `spice_net_list.sp` into ADS for simulation

## Input Parameters
- **bigY** - Admittance matrix (Np × Np × D) where Np = number of ports, D = frequency points
- **s** - Frequency data points in Laplace form (s = jω)

## Output
- **bigYfit** - Fitted admittance data (accuracy depends on iterations and model order)
- **rmserr** - Root mean square error of the fitted model
- **spice_net_list.sp** - Passivity-enforced SPICE netlist
- **Visualization plots** - Saved to `images/` directory

## State Space Model
The system follows the standard form:
```
          N   SERCijm
  Yij(s)=SUM(---------) +SERDij +s*SEREij ("ij" denotes element i,j)
         m=1 (s-SERAm)
```

## Testing
- Unit tests available in `test/` directory
- `test_rpdriver.py` - Tests pole-residue processing
- `test_vfdriver.py` - Tests vector fitting functionality

## Dependencies
- NumPy
- Matplotlib
- SciPy
- Python 3.x (tested up to Python 3.7)

## References
The implementation is based on research papers and algorithms from:
- Bjorn Gustavsen and Phil Reinhold's vector fitting work
- Various IEEE and academic publications on passivity enforcement

## Author
Jared Reichle
- Work email: reic2736@vandals.uidaho.edu
- Personal email: reichle.jared@gmail.com

## Created
May 1st, 2020 