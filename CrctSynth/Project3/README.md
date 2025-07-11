# Project 3: Passivity Enforcement

## Overview
This project implements passivity enforcement algorithms for state-space models of arbitrary N-port networks. It ensures that synthesized circuits are stable and passive for time-domain simulation.

## Description
The passivity enforcement algorithm processes vector-fitted state-space models and ensures they meet passivity criteria. This is crucial for time-domain stability in circuit simulation. The program uses Y-parameter (Admittance matrix) representation and includes eigenvalue analysis to verify passivity.

## Key Features
- **Passivity Enforcement**: Ensures circuit models are passive and stable
- **Eigenvalue Analysis**: Visualizes and analyzes system eigenvalues
- **SPICE Netlist Generation**: Creates simulation-ready netlists
- **Multi-port Support**: Handles arbitrary N-port networks
- **State-space Processing**: Works with A, C, D, E matrices from vector fitting

## Files
- `Project3.py` - Main passivity enforcement script
- `SERA.csv`, `SERC.csv`, `SERD.csv`, `SERE.csv` - State-space matrices
- `sPR.csv` - Frequency samples matrix
- `P3.sp` - Generated SPICE netlist
- `Eigenvalues.jpg` - Eigenvalue analysis plot
- `Passivity_Enforcement_Jared_R.pptx` - Presentation slides

## Usage
1. Ensure all required CSV files are in the working directory
2. Run `Project3.py`
3. The program will display passivity enforcement results
4. A SPICE netlist (`P3.sp`) will be generated
5. Import the netlist into ADS for time-domain testing

## Required Input Files
- `SERA.csv` - A portion of the state space model
- `SERC.csv` - C portion of the state space model  
- `SERD.csv` - D portion of the state space model
- `SERE.csv` - E portion of the state space model
- `sPR.csv` - Frequency samples matrix

## State Space Model
The system is represented as:
```
          N   SERCijm
  Yij(s)=SUM(---------) +SERDij +s*SEREij ("ij" denotes element i,j)
         m=1 (s-SERAm)
```

## Output
- **Eigenvalue Plot**: Shows system stability characteristics
- **SPICE Netlist**: Passivity-enforced circuit for simulation
- **Passivity Metrics**: Verification of enforced passivity
- **Visualization**: Plots saved to `images/` directory

## Dependencies
- NumPy
- Matplotlib
- SciPy
- Python 3.x (tested up to Python 3.7)

## Testing
- Time-domain stability testing can be performed in ADS
- The program provides passivity verification metrics
- Eigenvalue analysis confirms stability characteristics

## Author
Jared Reichle
- Work email: reic2736@vandals.uidaho.edu
- Personal email: reichle.jared@gmail.com

## Created
April 12th, 2020 