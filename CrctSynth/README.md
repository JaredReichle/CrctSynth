# ECE404: Advanced Circuit Synthesis and Vector Fitting

## Overview
This repository contains a comprehensive collection of projects focused on vector fitting algorithms, circuit synthesis, and passivity enforcement for high-frequency circuit design. These projects demonstrate the complete workflow from frequency domain data to simulation-ready SPICE netlists.

## Project Structure

### 📁 Project 1: Vector Fitting Algorithm
**Location**: `Project1/`
**Purpose**: Implements the core Vector Fitting algorithm for frequency domain data analysis

**Key Features**:
- Multi-port support (single-port and 4-port systems)
- 18th order system approximations
- Complex data processing capabilities
- Automatic pole-residue fitting

**Main Scripts**:
- `ex3.py` - Single port system example
- `ex4.py` - 4-port system example
- `Project1.py` - Main project orchestration
- `Vectorfitting.py` - Core algorithm implementation

**Usage**:
```bash
cd Project1
python ex3.py  # Single port example
python ex4.py  # 4-port example
```

### 📁 Project 2: Circuit Synthesis
**Location**: `Project2/`
**Purpose**: Converts pole-residue data into equivalent circuit components and SPICE netlists

**Key Features**:
- Pole-residue to circuit component conversion
- Automatic RL/RLRC branch selection
- SPICE netlist generation
- Multi-port network support

**Main Script**:
- `Project2.py` - Main synthesis script

**Usage**:
```bash
cd Project2
python Project2.py
```

### 📁 Project 3: Passivity Enforcement
**Location**: `Project3/`
**Purpose**: Ensures synthesized circuits are passive and stable for time-domain simulation

**Key Features**:
- Passivity enforcement algorithms
- Eigenvalue analysis and visualization
- State-space model processing
- Y-parameter network support

**Main Script**:
- `Project3.py` - Main passivity enforcement script

**Usage**:
```bash
cd Project3
python Project3.py
```

### 📁 Final Project: Advanced Vector Fitting with Passivity Enforcement
**Location**: `Final Project/`
**Purpose**: Complete workflow combining fast-relaxed vector fitting with passivity enforcement

**Key Features**:
- Fast-relaxed vector fitting algorithm
- Complete passivity enforcement workflow
- Real-time visualization and progress tracking
- Comprehensive testing suite
- S-parameter to SPICE netlist conversion

**Main Script**:
- `Final_Project.py` - Main orchestration script

**Usage**:
```bash
cd "Final Project"
python Final_Project.py
```

## Prerequisites

### Software Requirements
- **Python**: 3.x (tested up to Python 3.7)
- **NumPy**: For numerical computations
- **Matplotlib**: For plotting and visualization
- **SciPy**: For scientific computing functions

### Installation
```bash
pip install numpy matplotlib scipy
```

## Quick Start Guide

### 1. Vector Fitting (Project 1)
Start with the vector fitting examples to understand the core algorithm:
```bash
cd Project1
python ex3.py  # Single port system
```

### 2. Circuit Synthesis (Project 2)
Convert pole-residue data to circuit components:
```bash
cd Project2
python Project2.py
```

### 3. Passivity Enforcement (Project 3)
Ensure circuit stability:
```bash
cd Project3
python Project3.py
```

### 4. Complete Workflow (Final Project)
Run the complete pipeline:
```bash
cd "Final Project"
python Final_Project.py
```

## File Organization

### Images Directory
Each project includes an `images/` subdirectory where all generated plots and visualizations are automatically saved. This keeps the project directories clean and organized.

### Data Files
- **CSV files**: Contain pole-residue data and state-space matrices
- **MAT files**: MATLAB format data (e.g., `ex2_Y.mat`)
- **NPZ files**: NumPy compressed data files
- **SP files**: SPICE netlists for circuit simulation

## Output Files

### Generated Netlists
- `spice_net_list.sp` - Main SPICE netlist (Projects 2 and Final)
- `P3.sp` - Project 3 passivity-enforced netlist

### Visualization Files
All plots are saved to the respective `images/` directories:
- Frequency response plots
- Eigenvalue analysis
- Fitting accuracy comparisons
- Passivity verification plots

## Testing

### Unit Tests
The Final Project includes a comprehensive test suite:
```bash
cd "Final Project/test"
python test_rpdriver.py
python test_vfdriver.py
```

### Validation
- Compare fitted data with original frequency response
- Verify passivity through eigenvalue analysis
- Test SPICE netlists in circuit simulators (ADS, PSPICE)

## Common Issues and Solutions

### Python Version Compatibility
- Ensure Python 3.x is used (tested up to 3.7)
- Some scripts may require specific NumPy/SciPy versions

### File Path Issues
- Use relative paths within each project directory
- Ensure all required data files are present before running scripts

### Memory Issues
- Large frequency datasets may require significant memory
- Consider reducing model order for very large datasets

## Contributing

### Code Style
- Follow PEP 8 Python style guidelines
- Include docstrings for all functions
- Add comments for complex algorithms

### Testing
- Add unit tests for new functionality
- Validate results against known test cases
- Ensure backward compatibility

## Author and Contact

**Jared Reichle**
- Work email: reic2736@vandals.uidaho.edu
- Personal email: reichle.jared@gmail.com

## License
This project is part of ECE404 coursework at the University of Idaho. Please respect academic integrity and proper attribution when using or modifying these algorithms.

## References
- Gustavsen, B., & Semlyen, A. (1999). Rational approximation of frequency domain responses by vector fitting.
- Various IEEE publications on passivity enforcement and circuit synthesis
- Bjorn Gustavsen and Phil Reinhold's vector fitting implementations 