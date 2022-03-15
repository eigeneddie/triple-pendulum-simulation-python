# TRIPLE PENDULUM SIMULATION WITH PYTHON
This script is a demonstration for computational dynamics
from scratch, which means no specific libraries are used: just numpy, 
matplotlib.pyplot, a minor help from pandas to inspect matrices), and 
a direct implementation of multibody dynamics.

Understanding dynamics and its complexities can give you
a certain appreciation and awe to how software libraries with multibody 
dynamics physics engine actually work to generate various solid
mechanical system. This script, in the other hand, is merely hard
coded to the specified model on this README.md. 

System model looks a bit like what is shown in the image below. Both joints B and C have a torsional spring and rotational damper.

![Image of Triple Pendulum Model](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/triplePendulum1.png)

## 1. How to play with the code
 The scripts simulates a TRIPLE PENDULUM given the following user input
 1. Length of all 3 links [m]
 2. Mass of all 3 links [kg]
 3. Simulation time settings [s]
 4. External force elements (see README.md for system model)
 5. Initial conditions for each independent coordinates (in this case, 
    the angle of all links w.r.t. the vertical line) [rad]
 
 You may play with these variables on the USER INPUT PARAMETERS 
 section (lines just after importing modules).

 To run the script, simply enter the 
 following usage to the terminal: `python main_tp.py`
 
## 2. How the code works

The program consists of a main program `main_tp.py` and three sub-modules 
for different purposes.
1. `forceModule.py` contains functions for force elements (e.g. springs and dampers).
2. `constraintModuleTP.py` contains functions that computes constraint equations, jacobian matrices, and functions relating to these vectors/matrices such as position analysis and quadratic velocity terms.
3. `calcModulTP.py` contains useful mathematical utility functions such as rotational, converting local coordinates to the global coordinates, constructing mass matrix, rod inertia, etc. 

The calculation program is based on the augmented dynamics equation form. The program simply constructs the required matrix and vectors to find the second derivative of the system's generalized coordinates. A runge kutta integration algorithm is used to find the position for each generalized coordinates.

## 3. Example of results. 

Here is an example of the input parameters of the triple pendulum system.

```
link1, link2, link3          = 1, 1, 1 # [m]
mass1, mass2, mass3          = 1.2, 1.2, 2.0 #[kg]
timeStart, timeEnd, stepSize = 0.0, 10.0, 0.05 # [s] 
krB, crB, krC, crC           = 1.0, 1.0, 1.0, 1.0 # [Nm/rad], [Nms/rad]

# initial conditions
theta1Init, theta2Init, theta3Init = np.pi/4, np.pi/3, np.pi/2 #[rad]
```

With these parameters, we can evaluate certain variables of interest of the system such as the kinematics (joint velocity, position, etc.) and kinetics (joint reaction forces). 

Below are several examples for the variables of interest after running the script.

### a. Angular position of each link (theta)

Note: the angular position of each link are the independent coordinates. This is a 3-DOF system.

![theta](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/angular-position.png)


### b. Angular velocity of each link (omega)

![omega](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/angular-velocity.png)


### c. Angular acceleration of each link (alpha)

![alpha](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/angular-acceleration.png)


### d. Joint reaction forces on the Y direction (Force_y)

![force_y](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/join-reaction-force-y.png)


### e. X-coordinate position of link COG (Rx)

![Rx](https://github.com/eigeneddie/triple-pendulum-simulation-python/blob/main/img/X-coordinate-of-link-centers.png)
