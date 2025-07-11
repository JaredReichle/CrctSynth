# Project 1: Vector Fitting Algorithm

## Overview
This project implements the Vector Fitting algorithm in Python for frequency domain data analysis. It takes frequency response data and returns a transfer function that accurately represents the magnitude/phase relationship of the network.

## Description
The Vector Fitting Algorithm is a powerful tool for rational approximation of frequency domain data. This implementation handles 18th order systems with varying frequency data points (both complex and real) and finds the best fit for the data, delivering results that closely match the given poles-residues and input transfer function.

## Key Features
- **Multi-port Support**: Handles both single-port and 4-port systems
- **High-order Systems**: Supports 18th order approximations
- **Complex Data Processing**: Handles both real and complex frequency data
- **Automatic Fitting**: Finds optimal pole-residue representations
- **Visualization**: Generates plots showing fitting accuracy

## Files
- `Project1.py` - Main project script
- `Vectorfitting.py` - Core vector fitting implementation
- `Vector_fitting.py` - Alternative implementation
- `ex1.py`, `ex2.py` - Example scripts with predefined parameters
- `ex3.py` - Single port system example
- `ex4.py` - 4-port system example
- `MLtoPY.py` - MATLAB to Python conversion utilities
- `Test.py` - Comprehensive testing suite

## Usage
1. **For predefined examples**: Run `ex3.py` (single port) or `ex4.py` (4-port)
2. **For custom parameters**: Modify `ex1.py` or `ex2.py` with your desired parameters
3. **For testing**: Run `Test.py` for comprehensive validation

## Parameters
- **s** - Frequency samples: Random points representing system frequency response
- **f** - Frequency response f(s): Model system for vector fitting
- **d** - Offset: Starting magnitude value
- **h** - Slope: Fitted vector slope and system fitting characteristics

## Algorithm Components
- **Main Function**: Orchestrates the fitting process
- **vectfit_step**: Performs individual fitting iterations
- **calculate_residues**: Computes residue values
- **vectfit_auto**: Automated fitting with convergence checking

## Output
- Transfer function representation
- Pole-residue pairs
- Fitting accuracy metrics
- Visualization plots (saved to `images/` directory)

## Dependencies
- NumPy
- Matplotlib
- SciPy
- Python 3.x (tested up to Python 3.7)

## Author
Jared Reichle
- Work email: reic2736@vandals.uidaho.edu
- Personal email: reichle.jared@gmail.com

## Created
January 23rd, 2020 