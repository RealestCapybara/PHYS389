# PHYS389 Simulation Project

## Simulation and diagnostic plotting for a 3D connected n-pendulum system

I was intesrested in modelling chaotic systems, and a simple case for this is the connected double pendulum. I was interested in expanding this further and simulating any number of connected pendula in 3D, as a test of my abilities.

## Table of Contents

- [Usage](#usage)
- [Design](#design)
- [Structure](#structure)

## Usage

In order to simulate a system, run Simulator.py, and to plot a system run Plotter.py. The specific functionality is determined by config.toml. When Simulator.py is run, it produces a .pickle file with the filename chosen by a setting in config.toml, and when Plotter.py is run it loads a pickle file with the same filename function given by config.toml. Besides the filename, the settings for simulating are held in the system.sim table, and the settings for plotting are held in the system.plot table.

The system.sim table has three system-wide varables, g, the gravitational acceleration, timedelta, the time in seconds between iterations, and steps, the numer of iterations for the simulator to calculate. Besides these, there is an array of pendula tables. For each connected pendulum desired, the user should add another [[system.sim.pendula]] object (along with relevant attributes). In each table, the user should specify the length of the pendulum, its mass, its polar (theta) and azimuthal (phi) angles, its polar and azimuthal angular speed, and which coordinate basis it is in z-polar (true) or x-polar (false). See the [Design](#design) section for more details on how and why these values are used.

The system.plot table just contains boolean values and lists of boolean values. Each bool can be set to true or false to plot the associated value on a value over time plot. For example 'x' controls whether the cartesian positions of the pendula will be plotted, with the first bool controlling the x-values, the sec ond bool the y-values, and the third controlling the z-values.

## Design

The resulting system can simulate and plot any number of arbitrarily oriented connected pendulum objects - however for the sake of simplification the model does assume that all mass of a pendulum is centred at its end, and assumes a frictionless environment. These assumptions were necessary to be able to complete the project within the required time window.
