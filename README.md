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

$$\frac{\partial L}{\partial \mathbf{\dot{p}}_j} = \sum_{i=j}^N m_i \sum_{k=1}^i \mathbf{\dot{p}}_k = \sum_{i=j}^N m_i \sum_{k=1}^i \dot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \dot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k}$$



$$\frac{d}{dt} \frac{\partial L}{\partial \mathbf{\dot{p}}_j} = \sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} $$

$$ 0 = \sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} - g \hat{k}$$

As each partial derivative, this can be simplified by taking the inner product of the above with both $\partial_{\theta} \mathbf{p}$ and $\partial_{\phi} \mathbf{p}$, which results in the following,

$$ 0 = \sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} - g \hat{k} \cdot \frac{\partial \mathbf{p}_k}{\partial \theta_k} $$

$$ 0 = \sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} - g \hat{k} \cdot \frac{\partial \mathbf{p}_k}{\partial \phi_k} $$

Which can then be simplified into a matrix equation:

$$ 0 = A 
\begin{pmatrix}
\ddot{\theta}_1 \\
\ddot{\phi}_1 \\
\ddot{\theta}_2 \\
\ddot{\phi}_2 \\
. \\
. \\
. \\
\ddot{\theta}_N \\
\ddot{\phi}_N 
\end{pmatrix} + B 
\begin{pmatrix}
\dot{\theta}_1^2 \\
2 \dot{\theta}_1 \dot{\phi}_1 \\
\dot{\phi}_1^2 \\
\dot{\theta}_2^2 \\
2 \dot{\theta}_2 \dot{\phi}_2 \\
\dot{\phi}_2^2 \\
. \\
. \\
. \\
\dot{\theta}_N^2 \\
2 \dot{\theta}_N \dot{\phi}_N \\
\dot{\phi}_N^2
\end{pmatrix} +
\begin{pmatrix}
M_1 g \partial_{\theta 1} z_1\\
M_1 g \partial_{\phi 1} z_1\\
M_2 g \partial_{\theta 2} z_2\\
M_2 g \partial_{\phi 2} z_2\\
. \\
. \\
. \\
M_N g \partial_{\theta N} z_N\\
M_N g \partial_{\phi N} z_N
\end{pmatrix}$$

where $M_i$ is given by the sum $\sum_{j=i}^N m_j$ where $m_j$ are the masses of each pendulum. In the formula, A and B are block matrices given by the following formulae:

$$A = 
\begin{pmatrix}
M_1A_{11}        &M_2A_{12}        &.&.&.&M_{N-1}A_{1(N-1)}    &M_NA_{1N} \\
M_2A_{21}        &M_2A_{22}        &.&.&.&M_{N-1}A_{2(N-1)}    &M_NA_{2N} \\
.                &.                &.& & &.                    &. \\
.                &.                & &.& &.                    &. \\
.                &.                & & &.&.                    &. \\
M_{N-1}A_{(N-1)1}&M_{N-1}A_{(N-1)2}&.&.&.&M_{N-1}A_{(N-1)(N-1)}&M_NA_{(N-1)N}\\
M_N A_{N1}        &M_N A_{N2}      &.&.&.&M_N A_{(N-1)N}       &M_NA_{NN}\\
\end{pmatrix}$$

$$B = 
\begin{pmatrix}
M_1B_{11}        &M_2B_{12}        &.&.&.&M_{N-1}B_{1(N-1)}    &M_NB_{1N} \\
M_2B_{21}        &M_2B_{22}        &.&.&.&M_{N-1}B_{2(N-1)}    &M_NB_{2N} \\
.                &.                &.& & &.                    &. \\
.                &.                & &.& &.                    &. \\
.                &.                & & &.&.                    &. \\
M_{N-1}B_{(N-1)1}&M_{N-1}B_{(N-1)2}&.&.&.&M_{N-1}B_{(N-1)(N-1)}&M_NB_{(N-1)N}\\
M_N B_{N1}        &M_N B_{N2}      &.&.&.&M_N B_{(N-1)N}       &M_NB_{NN}\\
\end{pmatrix}$$

where,

$$A_{ij} = 
\begin{pmatrix}
\partial_{\theta j}\mathbf{p}_j\cdot\partial_{\theta i}\mathbf{p}_i & \partial_{\phi j}\mathbf{p}_j\cdot\partial_{\theta i}\mathbf{p}_i \\
\partial_{\theta j}\mathbf{p}_j\cdot\partial_{\phi i}\mathbf{p}_i & \partial_{\phi j}\mathbf{p}_j\cdot\partial_{\phi i}\mathbf{p}_i
\end{pmatrix}$$

$$B_{ij} = 
\begin{pmatrix}
\partial_{\theta j}^2\mathbf{p}_j\cdot\partial_{\theta i}\mathbf{p}_i & 2\partial_{\theta j}\partial_{\phi j}\mathbf{p}_j \cdot \partial_{\theta i}\mathbf{p}_i & \partial_{\phi j}^2\mathbf{p}_j\cdot\partial_{\theta i}\mathbf{p}_i \\
\partial_{\theta j}^2\mathbf{p}_j\cdot\partial_{\phi i}\mathbf{p}_i & 2\partial_{\theta j}\partial_{\phi j}\mathbf{p}_j \cdot \partial_{\phi i}\mathbf{p}_i & \partial_{\phi j}^2\mathbf{p}_j\cdot\partial_{\phi i}\mathbf{p}_i \\
\end{pmatrix}$$


