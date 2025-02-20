# Sharpe-Ratio-Portfolio-
Matlab code 
# MATLAB SHARPE RATIO PORTFOLIO

## Overview
This project implements an optimization algorithm in MATLAB to perform matrix calculations and optimize variables using gradient descent techniques. The primary goal is to minimize a cost function by iteratively adjusting the variables to achieve the optimal solution.

## Files
- `SharpeRatio.m`: The main script that executes the optimization process.
- `example_output`:  A folder containing the five figures generated by the script. These figures correspond to the key results presented in the article.  
   - Figure 1: Functional over iterations.  
   - Figure 2: Evolution of strategies (probabilities).  
   - Figure 3: Norm of constraints over iterations.  
   - Figure 4: Sharpe ratio over iterations.  
   - Figure 5: Efficient frontier.
     
## Requirements
- MATLAB (tested on version R2023a or higher).
- Optimization Toolbox for MATLAB (if required by certain functions).

## Setup
1. Download or clone this repository to your local machine.
2. Open the `SharpeRatio.m` script in MATLAB.
3. Ensure all files are in the same directory or add the appropriate paths.

## Usage
1. Open MATLAB.
2. Navigate to the directory containing the `SharpeRatio.m` script.
3. Run the `SharpeRatio.m` script to start the optimization process.
4. The script will generate output plots, including:
   - A functional progression plot showing the convergence over iterations.
   - Variable progression plots showing the changes in decision variables throughout the optimization process.

## Example Output
All generated figures are provided in the `example_output` folder as JPG files.
   - Figure 1: Functional over iterations.  
   - Figure 2: Evolution of strategies (probabilities).  
   - Figure 3: Norm of constraints over iterations.  
   - Figure 4: Sharpe ratio over iterations.  
   - Figure 5: Efficient frontier.

## Results
The optimization algorithm is expected to converge towards an optimal solution, the progress of the optimization can be tracked through the plots, which display how both the objective function and the variables evolve during each iteration.

## Authors
- Julio Bernardo Clempner Kerik
   - Contact: jclempnerk@ipn.mx
- Alin Andrei Carsteanu
   - Contact: acarsteanu@ipn.mx
- Lesly Lisset Ortiz Cerezo
   - Contact: lortizc1401@alumno.ipn.mx
