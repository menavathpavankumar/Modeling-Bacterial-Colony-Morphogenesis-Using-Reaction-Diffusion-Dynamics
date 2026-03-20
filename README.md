# Modeling-Bacterial-Colony-Morphogenesis-Using-Reaction-Diffusion-Dynamics

## Overview

This project simulates bacterial colony growth using a reaction–diffusion framework inspired by Ben-Jacob et al. (1994).
It models how nutrient transport, bacterial proliferation, and signaling interactions lead to complex spatial patterns.

## Model Description

The system consists of three coupled fields:

* Nutrient concentration c(x, y, t)
* Bacterial density n(x, y, t)
* Signaling molecule concentration s(x, y, t)

These are governed by nonlinear partial differential equations derived from mass balance and diffusion principles.

## Numerical Approach

* Spatial discretization using finite difference method
* Time integration using Forward Euler scheme
* Stability ensured using CFL condition
* Boundary conditions: Zero-flux (Neumann)

## Features

* Simulation of colony growth over time
* Nutrient depletion analysis
* Signal propagation modeling
* Parameter sensitivity analysis

## Results

* Radial colony expansion observed over time
* Nutrient depletion at colony center
* Heatmaps showing effect of diffusion parameters
* Growth curves showing diffusion-limited dynamics

## Tech Stack

* Python
* NumPy
* Matplotlib

## How to Run

```bash
pip install -r requirements.txt
python main.py
```

## Key Insights

* Bacterial motility significantly influences colony size and morphology
* Nutrient diffusion affects growth rate and spatial distribution
* Signal decay parameter has highest impact on system behavior

## Author

Menavath Pavan Kumar
