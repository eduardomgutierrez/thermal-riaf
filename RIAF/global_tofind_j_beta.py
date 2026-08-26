#!/usr/bin/env python

from global_eqs import *
from scipy.optimize import root_scalar, root

vOut = -lamda*np.sqrt(sqrdSoundVel(temp_i_Out, temp_e_Out))

logTiOut = np.log(temp_i_Out/iMMW)
logTeOut = np.log(temp_e_Out/eMMW)
logvOut = np.log(-vOut)

print("INITIAL PARAMETERS")
print()
print("bhMass = ", blackHoleMass)
print("AccRate = ", accRateOut/eddAccRate)
print("beta = ", beta)
print("alpha = ", alpha)
print("s = ", s)
print("delta = ", delta)
print("rOut = ", rOut)
print("vOut = ", vOut)
print("TiOut = ", temp_i_Out)
print("TeOut = ", temp_e_Out)
print("machOut = ", lamda)
print("angVelOut = ", alpha*sqrdSoundVel(temp_i_Out, temp_e_Out) /
      (-vOut*rOut*schwRadius)/keplAngVel(rOut*schwRadius))
print()

root = root_scalar(bounds_beta, bracket=[log10j0, log10j1], args=(
    logTiOut, logTeOut, logvOut), method='toms748', maxiter=30)
