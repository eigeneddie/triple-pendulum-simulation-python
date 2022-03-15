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

## How to play with the code
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
 



