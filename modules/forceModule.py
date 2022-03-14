import numpy as np
from modules.calcModuleTP import link2index as l2i

def torSpring(kr, qi, i, j, theta0):
    thetai = qi[l2i(i, "theta")]
    thetaj = qi[l2i(j, "theta")]
    deltaTheta = thetai-thetaj
    Q_SpringThetai = -kr*(deltaTheta-theta0)
    Q_SpringThetaj = kr*(deltaTheta-theta0)

    return Q_SpringThetai, Q_SpringThetaj

def torDamp(cr, qiDot, i, j):
    omegai = qiDot[l2i(i, "theta")]
    omegaj = qiDot[l2i(j, "theta")]
    deltaOmega = omegai-omegaj
    Q_DampThetai = -cr*deltaOmega
    Q_DampThetaj = cr*deltaOmega

    return Q_DampThetai, Q_DampThetaj