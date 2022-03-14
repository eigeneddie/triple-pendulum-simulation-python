# frequently used calculation functions
import numpy as np
import pandas as pd

def ATransformMatrix (theta): #A
    theta = float(theta)
    ATransform = np.array([[np.cos(theta), -np.sin(theta)], 
                            [np.sin(theta), np.cos(theta)]])
    return ATransform 

def ATransformMatrixTHETA (theta): #A_theta
    theta = float(theta)
    ATransformTHETA = np.array([[-np.sin(theta), -np.cos(theta)], 
                                [np.cos(theta), -np.sin(theta)]])
    return ATransformTHETA

def prettyMatVect(matVect):
    prettyMatVect = pd.DataFrame(matVect, columns =
                                ['R1x', 'R1y', 'theta1',
                                'R2x','R2y','theta2'])
    return prettyMatVect

def prettyMatVect2(matVect):
    prettyMatVect = pd.DataFrame(matVect, columns =
                                ['1', '2', '3', '4', '5',
                                '1', '2', '3', '4', '5'])
    return prettyMatVect

def local2global(qi, u_bar_iP, link_i):
    # To calculate Point of Interest positions in terms of global coordinates
    index_i = link2index(link_i, "x")
    Ri = np.array([qi[index_i], qi[index_i+1]]) 
    riP = Ri + np.matmul(ATransformMatrix(qi[index_i+2]), u_bar_iP)
    
    return riP

def inertiaRod (mass, length):
    Ic = 1/12*mass*length**2
    return Ic

def massMatrix(massVect):
    n = np.size(massVect)
    massMat = np.identity(n)*massVect
    return massMat

def link2index(link, string):
    if string == "x":
        index = 3*(link-1)
    elif string == "y":
        index = 3*(link-1)+1
    elif string == "theta":
        index = 3*(link-1)+2
    index = int(index)
    return index