'''
 Copyright (c) 2022, EBS/MakotoYu95
 Multidynamics program for educational purposes.
 
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
 following usage to the terminal: 
 
 python main_tp.py

 This script is a demonstration for computational dynamics
 from scratch, which means no specific libraries are used: just 
 pure implementation of multibody dynamics.

 Understanding dynamics and its complexities can give you
 a certain appreciation and awe to how software libraries with multibody 
 dynamics physics engine actually work to generate various solid
 mechanical system. This script, in the other hand, is merely hard
 coded to the specified model on the README.md.
 
 Happy hacking!!
'''

# import necessary modules
import numpy as np  
import matplotlib.pyplot  as plt
from modules import calcModuleTP       as calMod
from modules import constraintModuleTP as conMod
from modules import forceModule        as fMod
from modules.calcModuleTP import link2index as l2i, prettyMatVect2

# 1. === USER INPUT PARAMETERS (GLOBAL VARIABLES) ======
# ======================================================
# ======================================================
link1, link2, link3          = 1, 1, 1 # [m]
mass1, mass2, mass3          = 1.2, 1.2, 2.0 #[kg]
timeStart, timeEnd, stepSize = 0.0, 10.0, 0.05 # [s] 
krB, crB, krC, crC           = 1.0, 1.0, 1.0, 1.0 # [Nm/rad], [Nms/rad]

# initial conditions
theta1Init, theta2Init, theta3Init = np.pi/4, np.pi/3, np.pi/2 #[rad]
#=======================================================
#=======================================================

# POINTS OF INTEREST, LOCAL JOINTS
u_bar_1A  = np.array([ [0], [ link1/2] ])
u_bar_1B  = np.array([ [0], [-link1/2] ])
u_bar_2B  = np.array([ [0], [ link2/2] ])
u_bar_2C  = np.array([ [0], [-link1/2] ])
u_bar_3C  = np.array([ [0], [ link3/2] ])

# stopping criteria, gravity
epsilon = 0.000000000001
gravity = 9.81 # [m/s^2]

# Derived parameters
inertiaJ1 = calMod.inertiaRod(mass1, link1)
inertiaJ2 = calMod.inertiaRod(mass2, link2)
inertiaJ3 = calMod.inertiaRod(mass3, link3)
massVector = np.array([[mass1], [mass1], [inertiaJ1],
                       [mass2], [mass2], [inertiaJ2],
                       [mass3], [mass3], [inertiaJ3]])
mass_Matrix =  calMod.massMatrix(massVector) 
time = np.arange(timeStart, timeEnd, stepSize).T

# 2. DEFINE GENERALIZED COORDINATES - 9x1
n, nc  = 9, 6 # Generalized coordinates, number of Constraints
qi             = np.zeros((n,1))    # gen. position
qiDot          = np.zeros((n,1))    # gen. velocity
qiDotDot_lamda = np.zeros((n+nc,1)) # gen. acceleration
qi[l2i(1, "theta")] = theta1Init
qi[l2i(2, "theta")] = theta2Init
qi[l2i(3, "theta")] = theta3Init

# Memory
q_allTime       = np.zeros((np.size(time),  n))
v_allTime       = np.zeros((np.size(time),  n))
a_allTime       = np.zeros((np.size(time),  n))
FReact_allTime  = np.zeros((np.size(time), nc))
constraintVect  = np.zeros((nc,1))

#========== MAIN PROGRAM ==================
def mainProg():
    global qi, qiDot, qiDotDot_lamda
    timeNow  = timeStart
    iter = 0

    #3. NEWTON RHAPSON
    # KNOWN   : qi_indep, qiDot_indep (@ t)
    # FIND    : 1. qi_dep, qiDot_dep, all accelerations (@ t) 
    #           2. qi_indep, qiDot_indep (@ t+1)

    for timeID in range(np.size(time)):
        max_iteration = 50
        count = 0
        delta_qDep_norm = 1
        # a. Find dependent position
        while delta_qDep_norm > epsilon:

            Cq, Cq_dep, Cq_indep, constraintVect = config(qi)
            q_dep = np.concatenate((qi[0:2], qi[3:5], qi[6:8]), axis = 0) # position
            q_depNew, delta_qDep_norm = conMod.positionAnalysis(constraintVect,
                                                                Cq_dep, q_dep) 
            count = count + 1
            if (delta_qDep_norm<epsilon) or (count>max_iteration):
               break

       # b. Store q_dep in qi
        qi[0:2] = q_depNew[0:2]
        qi[3:5] = q_depNew[2:4]
        qi[6:8] = q_depNew[4:6]

        # c. Find dependent velocity
        qDot_indep  = np.concatenate((qiDot[2:3], qiDot[5:6], qiDot[8:9]), axis = 0) # velocity
        Cdi = np.dot(np.linalg.inv(-Cq_dep), Cq_indep)
        qDot_dep = np.dot(Cdi, qDot_indep)

        # d. Store qDot_dep in qiDot
        qiDot[0:2], qiDot[3:5], qiDot[6:8]= qDot_dep[0:2], qDot_dep[2:4], qDot_dep[4:6]
    
        # e. Find dependent acceleration (indep acc at the same time)
        qiDotDot_lamda = systemEquation(0, Cq, qi, qiDot)

        # f. Store everything
        q_allTime[timeID,:]      = qi.T
        v_allTime[timeID,:]      = qiDot.T
        a_allTime[timeID,:]      = qiDotDot_lamda[0:n].T
        FReact_allTime[timeID,:] = qiDotDot_lamda[n:n+nc].T

        # g. Calculate q, qdot, qdotdot independent @ t+1
        qi, qiDot = rungeKutta4_AtTimeNow( qi, qiDot, systemEquation, 
                                            stepSize, timeNow)
        iter = iter +1
        timeNow = timeNow + stepSize

    plt.figure(1)
    plt.plot(time, q_allTime[:, l2i(1, "theta")])
    plt.plot(time, q_allTime[:, l2i(2, "theta")])
    plt.plot(time, q_allTime[:, l2i(3, "theta")])
    plt.title('theta')
    plt.ylabel('angular position [rad]')
    plt.xlabel('time [s]')
    plt.grid(True)
    plt.legend(["theta 1", "theta 2", "theta 3"])

    plt.figure(2)
    plt.plot(time, v_allTime[:, l2i(1, "theta")])
    plt.plot(time, v_allTime[:, l2i(2, "theta")])
    plt.plot(time, v_allTime[:, l2i(3, "theta")])
    plt.title('omega')
    plt.ylabel('angular speed [rad/s]')
    plt.xlabel('time [s]')
    plt.grid(True)
    plt.legend(["omega 1", "omega 2", "omega 3"])
    
    plt.figure(3)
    plt.plot(time, q_allTime[:, l2i(1, "x")])
    plt.plot(time, q_allTime[:, l2i(2, "x")])
    plt.plot(time, q_allTime[:, l2i(3, "x")])
    plt.title('Rx (absolute horizontal position) of every link')
    plt.ylabel('position [m]')
    plt.xlabel('time [s]')
    plt.grid(True)
    plt.legend(["Rx 1", "Rx 2", "Rx 3"])
    
    plt.figure(4)
    plt.plot(time, -FReact_allTime[:, 1])
    plt.plot(time, -FReact_allTime[:, 3])
    plt.plot(time, -FReact_allTime[:, 5])
    plt.title('Joint reaction force y-axis [N]')
    plt.ylabel('Force [N]')
    plt.xlabel('time [s]')
    plt.grid(True)
    plt.legend(["Force 1", "Force 2", "Force 3"])

    plt.figure(5)
    plt.plot(time, a_allTime[:, l2i(1, "theta")])
    plt.plot(time, a_allTime[:, l2i(2, "theta")])
    plt.plot(time, a_allTime[:, l2i(3, "theta")])
    plt.title('Alpha')
    plt.ylabel('angular acceleration [rad/s/s]')
    plt.xlabel('time [s]')
    plt.grid(True)
    plt.legend(["alpha 1", "alpha 2", "alpha 3"])
    
    plt.show()
#==========================================

# IMPORTANT CALCULATION FUNCTIONS
def config(qi): #OKAY!
    r1A = calMod.local2global(qi, u_bar_1A, 1)
    r1B = calMod.local2global(qi, u_bar_1B, 1)
    r2B = calMod.local2global(qi, u_bar_2B, 2)
    r2C = calMod.local2global(qi, u_bar_2C, 2)
    r3C = calMod.local2global(qi, u_bar_3C, 3)

    # 4. CONSTRAINT EQUATION C
    constraintVect       = conMod.constraintEquation(r1A, r1B, r2B, r2C, r3C)

    # 5. JACOBIAN MATRIX Cq
    Cq, Cq_dep, Cq_indep = conMod.jacobianMatrix(qi, u_bar_1A, u_bar_1B,
                                                 u_bar_2B, u_bar_2C, u_bar_3C)

    return Cq, Cq_dep, Cq_indep, constraintVect

def systemEquation(t, Cq, qi, qiDot):
    # Construct MCq matrix (MASS MODULE)
    massSize = mass_Matrix.shape[0]
    constVSize = constraintVect.shape[0]
    matDim =  massSize + constVSize
    massAugmented = np.zeros((matDim, matDim))
    massAugmented[0:massSize, 0:massSize] = mass_Matrix
    massAugmented[massSize:matDim, 0:massSize] = Cq
    massAugmented[0:massSize, massSize:matDim] = np.transpose(Cq)

    # Construct QeQd vector (FORCE MODULE)
    Qe = np.zeros((massSize,1), dtype = float)
    
    # External Force from Weight
    Qe[l2i(1, "y")] = -mass1*gravity
    Qe[l2i(2, "y")] = -mass2*gravity
    Qe[l2i(3, "y")] = -mass3*gravity 
    
    # External Force from spring
    # -joint B (link 1&2)
    QSpring1B, QSpring2B = fMod.torSpring(krB, qi, 1, 2, 0)
    # -joint C (link 2&3)
    QSpring2C, QSpring3C = fMod.torSpring(krC, qi, 2, 3, 0)

    # External Force from damper
    # -joint B (link 1&2)
    QDamp1B, QDamp2B = fMod.torDamp(crB, qiDot, 1, 2)
    # -joint C (link 2&3)
    QDamp2C, QDamp3C = fMod.torDamp(crC, qiDot, 2, 3)
    
    Qe[l2i(1,"theta")]= QSpring1B + QDamp1B
    Qe[l2i(2,"theta")]= QSpring2B + QSpring2C + QDamp2B + QDamp2C
    Qe[l2i(3,"theta")]= QSpring3C + QDamp3C

    Qd1 = conMod.QdCalc1(qi, qiDot, u_bar_1A, 1)
    Qd2 = conMod.QdCalc2(qi, qiDot, u_bar_1B, u_bar_2B, 1, 2)
    Qd3 = conMod.QdCalc2(qi, qiDot, u_bar_2C, u_bar_3C, 2, 3)
    Qd = np.concatenate((-Qd1, Qd2, Qd3), axis = 0) #6x1
    
    QeAug = np.concatenate((Qe, Qd), axis = 0) #15x1
    mass_MatInverse = np.linalg.inv(massAugmented)
    qiDotDot_lamda = np.dot(mass_MatInverse, QeAug)
    return qiDotDot_lamda

def rungeKutta4_AtTimeNow(qi, qiDot, systemFunction, stepSize, timeNow):
    # This function works with ANY number of DOF
    x = np.array([qi[l2i(1, "theta")], 
                  qi[l2i(2, "theta")],
                  qi[l2i(3, "theta")]])

    xDot = np.array([qiDot[l2i(1, "theta")], 
                     qiDot[l2i(2, "theta")],
                     qiDot[l2i(3, "theta")]])

    y = np.concatenate((x, xDot), axis = 0)
    numberOfDOF = int(np.size(y)/2)
    
    # RungeKutta4
    t1  = timeNow
    Cq, _, _, _ = config(qi)
    f_1 = systemFunction(t1, Cq, qi, qiDot)
    k1  = np.zeros((np.size(y), 1))        
    for x    in range(numberOfDOF):
        k1[x]   = y[x+numberOfDOF]
        k1[x+numberOfDOF] = f_1[l2i(x+1, "theta")]
    
    t2  = t1+0.5*stepSize
    y2  = y + 0.5*k1*stepSize
    for i in range(numberOfDOF):
        qi   [l2i(i+1, "theta")] = y2[i]
        qiDot[l2i(i+1, "theta")] = y2[i+numberOfDOF]
    Cq, _, _, _ = config(qi)
    f_2 = systemFunction(t2, Cq, qi, qiDot)
    k2  = np.zeros((np.size(y), 1))           
    for x    in range(numberOfDOF):
        k2[x]   = y2[x+numberOfDOF]
        k2[x+numberOfDOF] = f_2[l2i(x+1, "theta")]
    
    t3  = t1+0.5*stepSize
    y3  = y + 0.5*k2*stepSize
    for i in range(numberOfDOF):
        qi   [l2i(i+1, "theta")] = y3[i]
        qiDot[l2i(i+1, "theta")] = y3[i+numberOfDOF]
    Cq, _, _, _ = config(qi)
    f_3 = systemFunction(t3, Cq, qi, qiDot)
    k3  = np.zeros((np.size(y), 1))           
    for x    in range(numberOfDOF):
        k3[x]   = y3[x+numberOfDOF]
        k3[x+numberOfDOF] = f_3[l2i(x+1, "theta")]
    
    t4  = t1+stepSize
    y4  = y + k3*stepSize
    for i in range(numberOfDOF):
        qi   [l2i(i+1, "theta")] = y4[i]
        qiDot[l2i(i+1, "theta")] = y4[i+numberOfDOF]
    Cq, _, _, _ = config(qi)
    f_4 = systemFunction(t4, Cq, qi, qiDot)
    k4  = np.zeros((np.size(y), 1))
    for x    in range(numberOfDOF):
        k4[x]   = y4[x+numberOfDOF]
        k4[x+numberOfDOF] = f_4[l2i(x+1, "theta")]

    RKFunct = (k1 + 2*k2 + 2*k3 + k4)/6

    yNew = y + stepSize*RKFunct
    for i in range(numberOfDOF):
        qi   [l2i(i+1, "theta")] = yNew[i]
        qiDot[l2i(i+1, "theta")] = yNew[i+numberOfDOF]

    return qi, qiDot

# Run main program
if __name__=="__main__":
    mainProg()