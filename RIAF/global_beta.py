#!/usr/bin/env python

from global_eqs import *

vOut = -lamda*np.sqrt(sqrdSoundVel(temp_i_Out,temp_e_Out))

logTiOut = np.log(temp_i_Out/iMMW)
logTeOut = np.log(temp_e_Out/eMMW)
logvOut = np.log(-vOut)

print("INITIAL PARAMETERS")
print()
print("bhMass = ",blackHoleMass)
print("AccRate = ",accRateOut/eddAccRate)
print("beta = ",beta)
print("alpha = ",alpha)
print("s = ",s)
print("delta = ",delta)
print("rOut = ",rOut)
print("vOut = ",vOut)
print("TiOut = ",temp_i_Out)
print("TeOut = ",temp_e_Out)
print("machOut = ",lamda)
print("angVelOut = ",alpha*sqrdSoundVel(temp_i_Out,temp_e_Out)/(-vOut*rOut*schwRadius)/keplAngVel(rOut*schwRadius))
print()

j0 = float(input("l_in = "))
y0 = np.array([logTiOut,logTeOut,logvOut,j0])

solution = solve_ivp(rhs_beta,(np.log(rOut),0.23),y0,method='LSODA',events=event_beta,dense_output=True)
print(solution.message)
print()
print("rSonic = ",np.exp(solution.t[-1]))
print()

ant = int(input("how many steps before the last one? "))

dlogr_sonic = solution.t[-1]-solution.t[-ant]
dlogTidlogr_sonic = (solution.y[0][-1]-solution.y[0][-ant])/(solution.t[-1]-solution.t[-ant])
logTi_sonic = solution.y[0][-1]+dlogTidlogr_sonic*dlogr_sonic
dlogTedlogr_sonic = (solution.y[1][-1]-solution.y[1][-ant])/(solution.t[-1]-solution.t[-ant])
logTe_sonic = solution.y[1][-1]+dlogTedlogr_sonic*dlogr_sonic
dlogvdlogr_sonic = (solution.y[2][-1]-solution.y[2][-ant])/(solution.t[-1]-solution.t[-ant])
logv_sonic = solution.y[2][-1]+dlogvdlogr_sonic*dlogr_sonic

logr_sonic = solution.t[-1]+dlogr_sonic
y0 = np.array([logTi_sonic,logTe_sonic,logv_sonic,j0])
solution2 = solve_ivp(rhs_beta,(logr_sonic,np.log(1.1)),y0,method='LSODA',vectorized=False,dense_output=True)
print(solution2.message)

logr = np.flip(np.concatenate((solution.t,solution2.t)))
logTi = np.flip(np.concatenate((solution.y[0],solution2.y[0])))
logTe = np.flip(np.concatenate((solution.y[1],solution2.y[1])))
logv = np.flip(np.concatenate((solution.y[2],solution2.y[2])))

r = np.exp(logr)*schwRadius
TiR = np.exp(logTi)*iMMW
TeR = np.exp(logTe)*eMMW
v = -np.exp(logv)
rho = massDensity(r,TiR,TeR,v)
cs = np.sqrt(sqrdSoundVel(TiR,TeR))
h = height(r,TiR,TeR)
logr10 = logr/np.log(10)
l = j0 - alpha*np.exp(logr)*(cs*cs)/cLight/v

gamma_i = 5.0/3.0
factor_gamma = np.sqrt( ( (3.0*gamma_i-1.0)+2.0*(gamma_i-1.0)*alpha*alpha )/(gamma_i+1.0) )

import matplotlib.pyplot as plt

plt.xlim([0,logr10[-1]])
plt.plot(logr10,np.log10(TiR),c='b')
plt.plot(logr10,np.log10(TeR),c='r')
plt.savefig("Temperatures.pdf")

plt.clf()
plt.xlim([0,logr10[-1]])
plt.ylim([0,2.2])
plt.plot(logr10,-v/cs*factor_gamma)
plt.savefig("MachNumber.pdf")

plt.clf()
plt.xlim([0,logr10[-1]])
plt.ylim([-3.5,-2])
plt.plot(logr10,np.log10(rho*2*h))
plt.savefig("SurfaceDens.pdf")

rad = np.logspace(0.01,4,1000) * schwRadius
plt.clf()
plt.xlim([0,logr10[-1]])
plt.ylim([0,1])
plt.plot(np.log10(rad/schwRadius),accRateRIAF(rad)/accRateOut, color='k')
plt.plot(np.log10(rad/schwRadius), gAux(rad), color='r')
plt.plot(np.log10(rad/schwRadius), fAux(rad), color='b')
plt.savefig("accRate.pdf")

plt.clf()
plt.xlim([0,logr10[-1]])
plt.ylim([0,1])
plt.plot(logr10,h/np.exp(logr)/schwRadius)
plt.savefig("HR.pdf")

plt.clf()
plt.xlim([0,logr10[-1]])
#plt.ylim()
plt.plot(logr10,np.log10(l))
plt.savefig("angularMom.pdf")

plt.clf()
plt.xlim([0,logr10[-1]])
plt.plot(logr10,np.log10(eDens(r,TiR,TeR,v)))
plt.savefig("eDens.pdf")

magf = np.sqrt(8.0*np.pi*(1.0-beta)*rho*cs*cs)
plt.clf()
plt.xlim([0,logr10[-1]])
plt.plot(logr10,np.log10(magf))
plt.savefig("magf.pdf")

np.savetxt(adafParameters,(blackHoleMass,accRateNorm,s,beta,alpha,j0,delta,innerRadiusSSD,pIndex))
np.savetxt(adafFile,np.column_stack([logr,logTi,logTe,logv]),fmt='%7.5f',header=str(np.size(logr)),comments='')
