# Project 2: Circuit Synthesis

## Overview
This project implements circuit synthesis from pole-residue data using the Vector Fitting algorithm. It converts frequency domain data into equivalent circuit components (resistors, inductors, capacitors) and generates SPICE netlists for simulation.

## Description
The program takes pole-residue data from CSV files and synthesizes equivalent circuit components. It handles both real and complex conjugate pole-residue pairs, automatically determining whether to use RL or RLRC branches based on the data characteristics.

## Key Features
- **Pole-Residue Analysis**: Processes complex frequency domain data
- **Circuit Synthesis**: Converts mathematical models to physical circuit components
- **SPICE Netlist Generation**: Creates simulation-ready netlists
- **Multi-port Support**: Handles N-port networks
- **Automatic Branch Selection**: Chooses appropriate RL or RLRC configurations

## Files
- `Project2.py` - Main synthesis script
- `fit_s.csv` - Input pole-residue data
- `fit_s_6.csv` - Additional test data
- `spice_net_list.sp` - Generated SPICE netlist
- `Circuit_synthesis_Jared_R.pptx` - Presentation slides

## Usage
1. Ensure your pole-residue data is in the correct CSV format
2. Run `Project2.py`
3. The script will generate a SPICE netlist file
4. Import the netlist into your circuit simulator (e.g., ADS, PSPICE)

## Input Format
The CSV file should contain:
- Row 0: Poles (complex numbers)
- Row 1+: Residues for each port (complex numbers)

## Output
- SPICE netlist with equivalent circuit components
- Each port synthesized as a subcircuit
- Components include resistors, inductors, and capacitors

## Dependencies
- NumPy
- Python 3.x

## Author
Jared Reichle
- Work email: reic2736@vandals.uidaho.edu
- Personal email: reichle.jared@gmail.com 