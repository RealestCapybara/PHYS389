# PHYS389 Simulation Project

## 3D N-Pendulum Simulator

This is software that can take an arbitrary connected chain of pendula in 3D, and simulate their movement over time, and also make diagnostic plots of the resulting simulation. I was intesrested in modelling chaotic systems, and a simple case for this is the connected double pendulum. I was interested in expanding this further and simulating any number of connected pendula in 3D, as a test of my abilities.

## Table of Contents

- [Usage](#usage)
- [Mathematics and Design](#mathematics-and-design)
- [Structure](#structure)
- [Future Work](#future-work)
- [Licence](#licence)

## Usage

In order to simulate a system, run `Simulator.py`, and to plot a system run `Plotter.py`. The specific functionality is determined by `config.toml`. When `Simulator.py` is run, it produces a .pickle file with the filename chosen by a setting in `config.toml`, and when `Plotter.py` is run it loads a pickle file with the same filename function given by `config.toml`. Besides the filename, the settings for simulating are held in the `system.sim` table, and the settings for plotting are held in the `system.plot` table.

The `system.sim` table has three system-wide varables, `g`, the gravitational acceleration, `timedelta`, the time in seconds between iterations, and `steps`, the numer of iterations for the simulator to calculate. Besides these, there is an array of pendula tables. For each connected pendulum desired, the user should add another `[[system.sim.pendula]]` item (along with relevant attributes). In each table, the user should specify the `length` of the pendulum, its `mass`, its polar (`theta`) and azimuthal (`phi`) angles, its polar and azimuthal angular speed, and which coordinate basis it is in z-polar (`true`) or x-polar (`false`). See the [Design](#design) section for more details on how and why these values are used.

The `system.plot` table just contains boolean values and lists of boolean values. Each bool can be set to `true` or `false` to plot the associated value on a value over time plot. For example `x` controls whether the cartesian positions of the pendula will be plotted, with the first bool controlling the x-values, the sec ond bool the y-values, and the third controlling the z-values. While most variables are variables associated with each pendulum, there are a few variables, namely `totalLM` (total linear momentum), `totalAM` (total angular momentum), `absTotalLM` (absolute total linear momentum), `absTotalAM` (absolute total angular momentum), `totalKE` (total kinetic energy), `totalPE` (total potential energy), and `totalME` (total mechanical energy), that are variables of the total system.

## Mathematics and Design

All of the following derivation was developed with reference to [Double pendulum](https://www.math.kyoto-u.ac.jp/~inou/preprint/doublependulum.pdf) by Hiroyuki Inou, 2018.

The resulting system can simulate and plot any number of arbitrarily oriented connected pendulum objects - however for the sake of simplification the model does assume that all mass of a pendulum is centred at its end, and assumes a frictionless environment. These assumptions were necessary to be able to complete the project within the required time window. It was also assumed that gravity is uniform and constant, and points in the negative z direction.

The differential equations used to construct the update functions used lagrangian mechanics, as it is preferable to newtonian mechanics when dealing with an arbitrary variable phase space. For each pendulum, only the polar angle, the azimuthal angle, and their time derivatives are necessary for constructing differential equations that can describe them.

In order to derive these differential equations, each pendulum was considered as an independent cartesian vector with x, y, and z all considered functions of the length of the pendulum, and the polar and azimuthal angles. However, as the length is constant in the case of an ideal pendulum, only the polar and azimuthal angles need to be considered - which helpfully reduces the number of independent variables necessary to consider. An issue with this construction, is that the points at which the polar angle (\$theta$) is 0 or pi, the angular speed of the azimuthal angle goes to infinity, which causes massive computational instability. A solution for this, is to construct a system whereby there are multiple spherical coordinate systems to be swapped between when sin(theta) is too small. This suggests two different methodologies, one whereby the entire connected pendulum system swaps between coordinate systems when one pendulum gets too close, and another whereby each pendulum swaps coordinate systems independently, and the system accounts for the possibility of interacting pendula in different spherical coordinate bases. The first possibility does limit the number of potential pendula that could be simulated, so my implimentation opts for the second methodology. Therefore in the derivation of the differential equations, the formulae of x, y, and z in terms of $\theta$ and $\phi$ are left unspecified. In the specific formulation chosen, the two spherical coordinate bases are a system with a z pole and one with an x pole. For convenience, the z-pole system is assumed as the default in most cases, and is defined such that $\theta \leq \pi/2$ gives you a negative z value.

Considering each pendulum as an independent vector then, the lagrangian can be constructed as:

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

$$ 0 = \left(\sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k} + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} - g \hat{k}\right) \cdot \frac{\partial \mathbf{p}_i}{\partial \theta_i} $$

$$ 0 = \left(\sum_{i=j}^N m_i \sum_{k=1}^i \ddot{\theta}_k \frac{\partial \mathbf{p}_k}{\partial \theta_k}  + \ddot{\phi}_k \frac{\partial \mathbf{p}_k}{\partial \phi_k} + \dot{\theta}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k^2} + 2 \dot{\theta}_k \dot{\phi}_k \frac{\partial^2 \mathbf{p}_k}{\partial \theta_k \partial \phi_k} + \dot{\phi}_k^2 \frac{\partial^2 \mathbf{p}_k}{\partial \phi_k^2} - g \hat{k}\right) \cdot \frac{\partial \mathbf{p}_i}{\partial \phi_i} $$

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

This matrix equation can then be rearranged to produce,

$$\begin{pmatrix}
\ddot{\theta}_1 \\
\ddot{\phi}_1 \\
\ddot{\theta}_2 \\
\ddot{\phi}_2 \\
. \\
. \\
. \\
\ddot{\theta}_N \\
\ddot{\phi}_N 
\end{pmatrix} = -A^{-1} \left( B 
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
\end{pmatrix} \right)$$

This can then be used to calculate the next position using multi-variable runge-kutte methods. This requires relabelling the first order time derivative of theta and phi as independent variables to produce the following equations,

$$\begin{pmatrix}
\dot{\omega}_{\theta 1} \\
\dot{\omega}_{\phi 1} \\
\dot{\omega}_{\theta 2} \\
\dot{\omega}_{\phi 2} \\
. \\
. \\
\dot{\omega}_{\theta N} \\
\dot{\omega}_{\phi N} 
\end{pmatrix} = -A^{-1} \left( B 
\begin{pmatrix}
\omega_{\theta 1}^2 \\
2 \omega_{\theta 1} \omega_{\phi 1} \\
\omega_{\phi 1}^2 \\
\omega_{\theta 2}^2 \\
2 \omega_{\theta 2} \omega_{\phi 2} \\
\omega_{\phi 2}^2 \\
. \\
. \\
. \\
\omega_{\theta N}^2 \\
2 \omega_{\theta N} \omega_{\phi N} \\
\omega_{\phi N}^2 \\
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
\end{pmatrix} \right)$$

and,

$$\begin{pmatrix}
\dot{\theta}_1 \\
\dot{\phi}_1 \\
\dot{\theta}_2 \\
\dot{\phi}_2 \\
. \\
. \\
\dot{\theta}_N \\
\dot{\phi}_N 
\end{pmatrix} = 
\begin{pmatrix}
\omega_{\theta 1} \\
\omega_{\phi 1} \\
\omega_{\theta 2} \\
\omega_{\phi 2} \\
. \\
. \\
\omega_{\theta N} \\
\omega_{\phi N} 
\end{pmatrix}$$

if these two equations are taken together as $\dot{y} = \text{f}(y)$, then RK4 can be applied.

## Structure

The core object used is the `Pendulum` object, which is specified in `Pendulum.py`. It has attributes for specifying `length`, `mass`, angular position (`pos`), angular velocity (`vel`) and `zPolar`. `zPolar` specifies whether the coordinate system has a z-polar coordinate system (`True`) or an x-polar coordinate system (`False`). Both angular position and velocity are 2 valued numpy arrays, where the first position is polar and the second is azimuthal. The Pendulum object specifies several partial derivative methods, which return the numerical partial derivatives of the pendulum position vector with respect to the polar angle, azimuthal angle, or some combination thereof. There are also methods for returning the numerical value of the cartesian position and velocity. Additionally, there is a method for fixing the coordinate basis, such that if the polar angle ($\theta$) is too close to a pole, the values of $\theta$, $\phi$, $\dot{\theta}$, and $\dot{\phi}$ are changed such that in the new coordinate basis the cartesian position and velocity is the same before and after, and also the value of zPolar is changed to indicate this. Additionally, this method that fixes the coordinates also ensures that $\theta$ and $\phi$ stay within the respective ranges $[0, \pi)$ and $[0, 2\pi)$, regardless of whether the coordinate basis is changed.

The secondary object used is the `Chain` object, which is specified in `Chain.py`. It has attributes specifying the `Pendulum` objects in the chain, and the value of the gravitational acceleration in the system. It has methods for calculating the next state of the system and for fixing the coordinates of the pendulum vectors. There are two methods for calculating the A block element and B block element from the inner products of the partial derivatives of the pendulum vectors with respect to theta, phi or some combination of the two. An additional method calculates the A and B matrices from these block elements. A further method calculates the Vector in the differential matrix equation that accounts for the gravitational acceleration. These class methods are used in a formulation of the RK4 method. See [Design](#design) for more details on the reason for these methods. In addition to this, there are class methods for calculating the cartesian position, velocity, linear momentum, angular momentum, and the kinetic, potential and mechanical energies of each pendulum - which are mostly used for plotting purposes.

These Objects are then used by Simulator.py and Plotter.py to simulate and plot the system. Both reference config.toml for the exact settings used. With an initial state, and filename, timedelta and number of steps, a system can be simulated for an arbitrary period of time. A smaller timedelta is more likely to avoid innacuracy, however, a smaller timedelta requires more steps for the same time period, and so will generally require more intensive calculation. The attributes of the pendulum objects are specified with an array of tables in the config.toml file. For each tick of the simulation, the chain object at that moment is saved to a list, which is then pickled, along with the timestep. Then during plotting, this list is extracted and processed into all of the potential data streams that are desired to be plotted in the config.toml file. While an animation would be very much possible, as it is not necessary for simulation or testing purposes, the Plotter.py file only has capability to make 2D diagnostic plots of the variable's value over time.

Lastly, `testing.py` is a simple unit testing file that runs checks to ensure that `Chain` and `Pendulum` methods return their expected values. It runs without failing, which means that the Chain and Pendulum objects are built to the design goals. If this project expands to other areas, then this testing file should be expanded to account for that.

## Future Work

While I am relatively satisfied with the program in its current state, there are a few sub-areas that it feels that this program could be expanded to encompass. Firstly, I would like to rework the system such that it could account for mass not centred at the end of the pendulum - initial ideas are a floating point value between 0 and 1 that determines the centre of mass of the pendulum, where 1 places the mass at the end of the pendulum vector. 

Additionally, I would like for there to be a way to account for the possibility of damping on the spherical pendulum joints, and general friction due to viscosity of the medium the system is placed within. I would also like for the system to be able to account for imperfect spherical joints, perhaps with more resistance in some directions, or maybe fully restricted to one particular axis. Further, extending the system to also account for connected spring systems would be interesting, as spring-pendulum systems give some interesting effects. This would require an implementation of non-constant length and hooke's law. Given that the length and angular position are disconnected variables in the Pendulum object, this feels like a possible extension - but it would of course require thorough investigation. 

Lastly, I would like for more complex forces on the system to be accounted for, perhaps a more interesting, or potentially non-constant force field to account for more environmental effects like wind or external movement. This extension of the forces could even potentially implement electromagnetic effects, like the lorentz force; or even implement diamagnetic, paramagnetic, and ferromagnetic materials in the makeup of the pendulum system.

## Licence

### The MIT License (MIT)

Copyright © 2023 Owen Wray

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
