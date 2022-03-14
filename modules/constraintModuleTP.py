#equations for constraints
import numpy        as np
from   modules.calcModuleTP   import ATransformMatrixTHETA as A_Theta, link2index
from   modules.calcModuleTP   import ATransformMatrix      as A_i

def constraintEquation(r1A, r1B, r2B, r2C, r3C):
    
    constraintVector = np.zeros((6,1))

    # Pin joint A
    constraintPinA = -r1A 
    for i in range(np.size(constraintPinA)):
        # Equation 1-2
        constraintVector[i] = constraintPinA[i]

    # Pin joint B
    constraintPinB = revolutJoint(r1B, r2B)
    for i in range(np.size(constraintPinB)):
        # Equation 3-4
        constraintVector[i+2] = constraintPinB[i]

    # Pin joint C
    constraintPinC = revolutJoint(r2C, r3C)
    for i in range(np.size(constraintPinC)):
        # Equation 3-4
        constraintVector[i+4] = constraintPinC[i]

    return constraintVector

def jacobianMatrix(qi, u_bar_1A, u_bar_1B, u_bar_2B, u_bar_2C, u_bar_3C):
    
    genCoor = np.size(qi) # number of generalized coordinates
    constEq = 6 # number of constraint equations

    jacobianMatrixCq = np.zeros((constEq, genCoor))
    identity2x2 = np.identity(2)

    # row 1-2
    Cq12 = np.dot(A_Theta(qi[link2index(1,"theta")]), u_bar_1A)
    jacobianMatrixCq[0:2,0:2] = -identity2x2
    jacobianMatrixCq[0:2,2:3] = -Cq12

    # row 3-4 (r1A = r2A)
    Cq34_link1 = np.dot(A_Theta(qi[link2index(1,"theta")]), u_bar_1B)
    Cq34_link2 = np.dot(A_Theta(qi[link2index(2,"theta")]), u_bar_2B)

    jacobianMatrixCq[2:4,0:2] = identity2x2
    jacobianMatrixCq[2:4,2:3] = Cq34_link1
    jacobianMatrixCq[2:4,3:5] = -identity2x2
    jacobianMatrixCq[2:4,5:6] = -Cq34_link2
    
    # row 5-6 (r2C = r3C)
    Cq56_link2 = np.dot(A_Theta(qi[link2index(2,"theta")]), u_bar_2C)
    Cq56_link3 = np.dot(A_Theta(qi[link2index(3,"theta")]), u_bar_3C)

    jacobianMatrixCq[4:6,3:5] = identity2x2
    jacobianMatrixCq[4:6,5:6] = Cq56_link2
    jacobianMatrixCq[4:6,6:8] = -identity2x2
    jacobianMatrixCq[4:6,8:9] = -Cq56_link3

    # SLICING
    # a. jacobian dependent
    jacobian_dependent = np.concatenate((jacobianMatrixCq[:,0:2], 
                                         jacobianMatrixCq[:,3:5],
                                         jacobianMatrixCq[:,6:8]), axis = 1)
    # b. jacobian independent
    jacobian_independent = np.concatenate((jacobianMatrixCq[:,2:3], 
                                           jacobianMatrixCq[:,5:6],
                                           jacobianMatrixCq[:,8:9]), axis = 1)

    return jacobianMatrixCq, jacobian_dependent, jacobian_independent

def positionAnalysis(constraintVector, jacobianMatrix, qi):
    inverse_jacobian = np.linalg.inv(jacobianMatrix)
    delta_qi = - np.matmul(inverse_jacobian, constraintVector)
    delta_qi_norm = np.linalg.norm(delta_qi)
    qi = qi + delta_qi

    return qi, delta_qi_norm

def QdCalc1(qi, qiDot, u_bar_iP, i):
    id = link2index(i, "theta")
    Qd = np.square(float(qiDot[id]))*np.dot(A_i(qi[id]), u_bar_iP)
    return Qd 

def QdCalc2(qi, qiDot, u_bar_iP, u_bar_jP, i, j): 
    id = link2index(i, "theta")
    jd = link2index(j, "theta")
    Qda = np.square(float(qiDot[id]))*np.dot(A_i(qi[id]), u_bar_iP)
    Qdb = np.square(float(qiDot[jd]))*np.dot(A_i(qi[jd]), u_bar_jP)
    Qd = Qda-Qdb
    return Qd

def revolutJoint (riP, riJ):
    constraintPin = riP-riJ
    
    return constraintPin