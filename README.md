# PHYS389 Simulation Project

## Simulation and diagnostic plotting for a 3D connected arbitrary n-pendulum system

I was intesrested in modelling chaotic systems, and a simple case for this is the connected double pendulum. I was interested in expanding this further and simulating any number of connected pendula in 3D, as a test of my abilities.

## Table of Contents

- [Usage](#usage)
- [Design](#design)
- [Structure](#structure)

## Usage

In order to simulate a system, run Simulator.py, and to plot a system run Plotter.py. The specific functionality is determined by config.toml. When Simulator.py is run, it produces a .pickle file with the filename chosen by a setting in config.toml, and when Plotter.py is run it loads a pickle file with the same filename function given by config.toml. Besides the filename, the settings for simulating are held in the system.sim table, and the settings for plotting are held in the system.plot table.

The system.sim table has three system-wide varables, g, the gravitational acceleration, timedelta, the time in seconds between iterations, and steps, the numer of iterations for the simulator to calculate. Besides these, there is an array of pendula tables. For each connected pendulum desired, the user should add another [[system.sim.pendula]] object (along with relevant attributes). In each table, the user should specify the length of the pendulum, its mass, its polar (theta) and azimuthal (phi) angles, its polar and azimuthal angular speed, and which coordinate basis it is in z-polar (true) or x-polar (false). See the [Design](#design) section for more details on how and why these values are used.

The system.plot table just contains boolean values and lists of boolean values. Each bool can be set to true or false to plot the associated value on a value over time plot. For example 'x' controls whether the cartesian positions of the pendula will be plotted, with the first bool controlling the x-values, the sec ond bool the y-values, and the third controlling the z-values. While most variables are variables associated with each pendulum, there are a few variables, namely totalLM (total linear momentum), totalAM (total angular momentum), absTotalLM (absolute total linear momentum), absTotalAM (absolute total angular momentum), totalKE (total kinetic energy), totalPE (total potential energy), and totalME (total mechanical energy), that are variables of the total system.

## Design

The resulting system can simulate and plot any number of arbitrarily oriented connected pendulum objects - however for the sake of simplification the model does assume that all mass of a pendulum is centred at its end, and assumes a frictionless environment. These assumptions were necessary to be able to complete the project within the required time window.

The differential equations used to construct the update functions used lagrangian mechanics, as it is preferable to newtonian mechanics when dealing with an arbitrary variable phase space. For each pendulum, only the polar angle, the azimuthal angle, and their time derivatives are necessary for constructing differential equations that can describe them.

In order to derive these differential equations, each pendulum was considered as an independent cartesian vector with x, y, and z all considered functions of the length of the pendulum, and the polar and azimuthal angles. However, as the length is constant in the case of an ideal pendulum, only the polar and azimuthal angles need to be considered - which helpfully reduces the number of independent variables necessary to consider. Considering each pendulum as an independent vector then, the lagrangian can be constructed as:

$$L = \frac{1}{2} \sum_{i=1}^{N} m_i \left( \sum_{j=1}^{i} \mathbf{\dot{p}}_j \right)^2 + g \sum_{i=1}^N m_i \left(\sum_{j=1}^i z_j\right)$$

Where $\mathbf{\dot{p}}_j$ are the time derivatives of the pendulum vectors, and $z_j$ are the cartesian z-components of the pendulum vectors, and $m_i$ are the masses of each pendulum. From this, the lagrangian equations can be derived:

$$\frac{\partial L}{\partial \mathbf{p}_j} = -g \sum_{i=j}^N 
\begin{pmatrix} 0 \\
0\\
m_i 
\end{pmatrix}$$

$$\frac{\partial L}{\partial \mathbf{\dot{p}}_j} = \sum_{i=j}^N m_i \sum_{k=1}^i \mathbf{\dot{p}}_k = \sum_{i=j}^N m_i \sum_{k=1}^i \theta_k \frac{\partial \mathbf{\dot{p}}_j}{\partial \theta_k} + \phi_k \frac{\partial \mathbf{\dot{p}}_j}{\partial \phi_k}$$


