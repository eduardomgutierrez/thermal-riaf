import scipy.special as sp
import scipy.integrate as integ

from scipy.integrate import solve_ivp
from astropy.modeling.models import BlackBody
import astropy.units as u
from constants_cgs import *

from global_initCond import *

# CONSTANTS [CGS]

mu = const.physical_constants['atomic mass constant'][0]*1.0e3
re = const.physical_constants['classical electron radius'][0]*1.0e2
eMMW = 1.14
iMMW = 1.23

# ADAF EQS

schwRadius = 2.0*G*blackHoleMass*solarMass / cLight2

def a_aux(x):
    k1 = sp.kn(1, x)
    k2 = sp.kn(2, x)
    k3 = sp.kn(3, x)
    result = np.where(k2 > 0.0, x*(0.25*(3.0*k3+k1)/k2-1.0), 1.5053)
    return result

def da_aux(x):
    k0 = sp.kn(0, x)
    k1 = sp.kn(1, x)
    k2 = sp.kn(2, x)
    k3 = sp.kn(3, x)
    k4 = sp.kn(4, x)
    aux = 0.125*(3.0*k4+4.0*k2+k0)
    result_aux = a_aux(x)/x - x*(aux - 0.5*(k1+k3)*(a_aux(x)/x+1.0))/k2
    return np.where(k2 > 0.0, result_aux, a_aux(x)/x)

def aux_der(x):
    return - da_aux(1.0/x) / (x*x)

def keplAngMom(rnorm):   # [rS c]
    return np.power(rnorm, 1.5)/(rnorm-1.0) / np.sqrt(2.0)

def keplAngVel(r):
    return np.sqrt(G*blackHoleMass*solarMass/r) / (r-schwRadius)

def dlogOmK_dlogr(rnorm):
    return -0.5 - rnorm/(rnorm-1.0)

def sqrdSoundVel(temp_i, temp_e):
    return boltzmann/(beta*mu) * (temp_i/iMMW + temp_e/eMMW)

def height(r, temp_i, temp_e):
    return np.sqrt(sqrdSoundVel(temp_i, temp_e))/keplAngVel(r)

def fAux(r):
    rTr = innerRadiusSSD * schwRadius
    rOutADAF = rOut * schwRadius
    return (1.0 - np.power(rTr/r,pIndex)) / (1.0 - np.power(rTr/rOutADAF,pIndex))

def gAux(r):
    rOutADAF = rOut * schwRadius
    rTr = innerRadiusSSD * schwRadius 
    return ((-(r**(pIndex + s) - rOutADAF**(pIndex + s)))*rTr**pIndex*s -
            rOutADAF**s*(r*rOutADAF)**pIndex*((rTr/r)**pIndex -
            (r/rOutADAF)**s*(rTr/rOutADAF)**pIndex)*(pIndex+s)) / (rOutADAF**s*(r*rOutADAF)**pIndex)/((-1 + (rTr/rOutADAF)**pIndex)*(pIndex + s))

def accRateRIAF(r):
    rTr = innerRadiusSSD * schwRadius
    rOutADAF = rOut * schwRadius
    return np.where(SSDdisk == 0, accRateOut * np.power(r/rOutADAF, s),
            np.where( r < rTr, accRateOut * gAux(rTr) * np.power(r/rTr, s), accRateOut * ( correctorAccRate * (1.0-rTr/r) + gAux(r) ) ) )


#    return np.where(SSDdisk == 0 or innerRadiusSSD < 5.0, accRateOut * np.power(r/rOutADAF, s),
#            np.where( r < rTr, accRateOut * gAux(rTr) * np.power(r/rTr, s), accRateOut * ( correctorAccRate * (1.0-rTr/r) + gAux(r) ) ) )

def massDensity(r, temp_i, temp_e, radialVel):
    return accRateRIAF(r) / (4.0*np.pi*r*height(r, temp_i, temp_e)*(-radialVel))

def magField(r, temp_i, temp_e, radialVel):
    return np.sqrt((1.0-beta)*massDensity(r, temp_i, temp_e, radialVel)*sqrdSoundVel(temp_i, temp_e)*8.0*np.pi)

def iDens(r, temp_i, temp_e, radialVel):
    return massDensity(r, temp_i, temp_e, radialVel)/(iMMW*mu)

def eDens(r, temp_i, temp_e, radialVel):
    return massDensity(r, temp_i, temp_e, radialVel)/(eMMW*mu)

def qie_exact(r, temp_i, temp_e, v):
    lnLambda = 20.0
    theta_e = boltzmann*temp_e/(electronMass*cLight2)
    theta_i = boltzmann*temp_i/(protonMass*cLight2)
    xe = 1.0/theta_e
    xi = 1.0/theta_i
    xei = xe+xi
    k2i = sp.kn(2, xi)
    k2e = sp.kn(2, xe)
    k1ei = sp.kn(1, xei)
    k0ei = sp.kn(0, xei)
    sumtheta = theta_i + theta_e
    aux1 = 1.875*thomson*(electronMass/protonMass)*cLight * lnLambda * eDens(r, temp_i, temp_e, v)*iDens(r, temp_i, temp_e, v) *\
        boltzmann*(temp_i-temp_e)
    aux2 = (2.0*sumtheta*sumtheta+1.0)/sumtheta
    return np.where(xei > 300.0,
                    np.where(xi > 150.0,
                             np.where(xe > 150.0,
                                      aux1*np.sqrt(2.0*xei /
                                                   (np.pi*xe*xi))*(aux2+2.0),
                                      aux1*np.sqrt(xi/xei)*(aux2+2.0)*np.exp(-xe)/k2e),
                             aux1*np.sqrt(xe/xei)*(aux2+2.0)*np.exp(-xi)/k2i),
                    aux1*(aux2 * k1ei/k2i + 2.0*k0ei/k2i)/k2e)

def qie_approx(r, temp_i, temp_e, v):
    theta_e = boltzmann*temp_e / (electronMass*cLight2)
    theta_i = boltzmann*temp_i / (protonMass*cLight2)
    lnLambda = 20.0
    diftemps = (boltzmann*temp_i - boltzmann*temp_e)
    aux = 1.5*electronMass/protonMass * eDens(r, temp_i, temp_e, v)*iDens(r,
                                                          temp_i, temp_e, v)*thomson*cLight*lnLambda*diftemps
    aux2 = (np.sqrt(2.0/np.pi)+np.sqrt(theta_e+theta_i)) / \
        np.power(theta_e+theta_i, 1.5)
    return aux*aux2

def gaunt(nu, temp_e):
    aux = boltzmann*temp_e/(planck*nu)
    zeda = 0.57695
    return np.where(aux < 1.0, np.sqrt(3.0/np.pi * aux), np.sqrt(3)/np.pi*np.log(4.0/zeda * aux))

def qei(r, temp_i, temp_e, v):
    theta = boltzmann*temp_e/(electronMass*cLight2)
    ne = eDens(r, temp_i, temp_e, v)
    ni = iDens(r, temp_i, temp_e, v)
    Fei = np.where(theta < 1.0,
                   4.0 * np.sqrt(2.0*theta/(np.pi*np.pi*np.pi)) *
                   (1.0+1.781*np.power(theta, 1.34)),
                   4.5*theta/np.pi * (np.log(1.123*theta+0.48)+1.5))
    return ne*ni*thomson*cLight*const.alpha*electronMass*cLight2*Fei

def qee(r, temp_i, temp_e, v):
    theta = boltzmann*temp_e/(electronMass*cLight2)
    ne = eDens(r, temp_i, temp_e, v)
    return ne*ne*cLight*re*re*const.alpha*electronMass*cLight2 * \
        np.where(theta < 1.0, 20.0/(9.0*np.sqrt(np.pi)) * (44.0-3.0*np.pi*np.pi) *
                 np.power(theta, 1.5)*(1.0+1.1*theta+theta*theta -
                                       1.25*np.power(theta, 2.5)),
                 24.0*theta*(np.log(1.1232*theta)+1.28))

def xiBremss(nu, r, temp_i, temp_e, v):
    qBremss = qee(r, temp_i, temp_e, v)+qei(r, temp_i, temp_e, v)
    return qBremss*np.exp(-planck*nu/(boltzmann*temp_e))*(planck/(boltzmann*temp_e))*gaunt(nu, temp_e)

def Iprim(x):
    return 4.0505/np.power(x, 1.0/6.0) * (1.0+0.4/np.power(x, 0.25)+0.5316/np.sqrt(x)) * \
        np.exp(-1.8899*np.power(x, 1.0/3.0))

def xiSync(nu, r, temp_i, temp_e, v):
    theta = boltzmann*temp_e/(electronMass*cLight2)
    nu0 = electronCharge*magField(r, temp_i, temp_e, v)/(2.0*np.pi*electronMass*cLight)
    xM = 2.0*nu/(3.0*nu0*theta*theta)
    k2 = sp.kn(2, 1.0/theta)
    k2 = np.where(k2 > 0, k2, 1.0e100)
    return 4.43e-30*4.0*np.pi*eDens(r, temp_i, temp_e, v)*nu / k2 * Iprim(xM)
    # return np.where(k2>0.0,4.43e-30*4.0*np.pi*eDens(r,temp_i,temp_e,v)*nu / sp.kn(2,1.0/theta) * \
    #        Iprim(xM),0.0)

# General functions

def bb_Intensity(freq,temp):
    return np.where(temp>0.0, 2*planck*np.power(freq,3)/cLight2 / (np.exp(planck*freq/(boltzmann*temp))-1.0), 0.0)

# SSD functions

def SSD_qPlus(radius):
    """Energy liberated per unit time per unit surface
    area from each face of the disk"""
    innerRadius = innerRadiusSSD*schwRadius
    const = 3*gravConstant*blackHoleMass*solarMass*accRateCD*fAux(radius) / (8*np.pi*np.power(radius,3))
    return np.where(radius > innerRadius, const * (1.0 - np.sqrt(innerRadius/radius)), 0.0)

def SSD_temp(radius):
    """Disk surface  temperature"""
    return np.power(SSD_qPlus(radius)/stefanBoltzmann,0.25)

def kappa(nu, r, temp_i, temp_e, v):
    bb = BlackBody(temp_e * u.K)
    bnu = 4.0*np.pi*bb(nu*u.Hz).value
    xi = xiBremss(nu, r, temp_i, temp_e, v)+xiSync(nu, r, temp_i, temp_e, v)
    return np.where(bnu > 1.0e-2*xi*height(r, temp_i, temp_e), xi/bnu, 50.0)

def tau(nu, r, temp_i, temp_e, v):
    return np.sqrt(np.pi)/2.0 * kappa(nu, r, temp_i, temp_e, v)*height(r, temp_i, temp_e)

def eta_esin(nu, r, tempi, tempe, v):

    theta = boltzmann*tempe/(electronMass*cLight2)
    ne = eDens(r, tempi, tempe, v)
    h = height(r, tempi, tempe)
    tau_es = 2.0*ne*thomson*h
    s_exp = tau_es*(tau_es+1.0)
    A = 1.0+4.0*theta*(1.0+4.0*theta)
    etaMax = 3.0*boltzmann*tempe/(planck*nu)
    jm = np.log(etaMax)/np.log(A)

    gamma1 = sp.gammainc(jm+1.0, A*s_exp)
    gamma2 = sp.gammainc(jm+1.0, s_exp)
    aux2 = s_exp*(A-1.0)

    result = np.where(aux2 < 200.0, np.exp(
        aux2)*(1.0-gamma1)+etaMax*gamma2, etaMax*gamma2)
    return result


def eta_dermer(nu, r, tempi, tempe, v):

    theta = boltzmann*tempe/(electronMass*cLight2)
    ne = eDens(r, tempi, tempe, v)
    h = height(r, tempi, tempe)
    tau_es = 2*ne*thomson*h
    P = 1.0-np.exp(-tau_es)
    x = planck*nu/(electronMass*cLight2)
    A = 1.0 + 4.0*theta*(1.0+4.0*theta)
    kap = P*(A-1.0)/(1.0-P*A)
    phi = -np.log(P)/np.log(A)
    x = x/(3.0*theta)
    result = 1.0 + kap*(1.0-np.power(x, phi-1.0))

    return np.where(result > 1.0, result, 1.0)

def flux(nu, r, temp_i, temp_e, v):
    bb = BlackBody(temp_e * u.K)
    bnu = bb(nu*u.Hz).value
    tauu = tau(nu, r, temp_i, temp_e, v)
    return 2.0*np.pi/np.sqrt(3.0) * bnu * (1.0-np.exp(-2.0*np.sqrt(3.0)*tauu))

def integrand_dermer(nu, r, temp_i, temp_e, v):
    result = flux(nu, r, temp_i, temp_e, v) * \
        eta_dermer(nu, r, temp_i, temp_e, v)
    result = np.where(result > 0.0, result, 0.0)
    return result


def integrand_esin(nu, r, temp_i, temp_e, v):
    result = flux(nu, r, temp_i, temp_e, v)*eta_esin(nu, r, temp_i, temp_e, v)
    result = np.where(result > 0.0, result, 0.0)
    return result

def qminus_esin(r, temp_i, temp_e, v):
    nuMin = 1.0e6
    nuMax = 1.0e21
    iter = 30
    paso = np.power(nuMax/nuMin, 1.0/float(iter-1))
    sum = 0.0
    nu = nuMin
    for i in range(iter):
        dnu = nu*(np.sqrt(paso)-1.0/np.sqrt(paso))
        sum += integrand_esin(nu, r, temp_i, temp_e, v)*dnu
        nu = nu*paso

    return sum / height(r, temp_i, temp_e)

def qminus_dermer(r, temp_i, temp_e, v):
    nuMin = 1.0e6
    nuMax = 1.0e21
    iter = 30
    paso = np.power(nuMax/nuMin, 1.0/float(iter-1))
    sum = 0.0
    nu = nuMin
    for i in range(iter):
        dnu = nu*(np.sqrt(paso)-1.0/np.sqrt(paso))
        sum += integrand_dermer(nu, r, temp_i, temp_e, v)*dnu
        nu = nu*paso
    return sum / height(r, temp_i, temp_e)

def qminus_dermer_fast(r, tempi, tempe, v):
    nuMin = 1.0e6
    nuMax = 1.0e21
    iter = 30
    paso = np.power(nuMax/nuMin, 1.0/float(iter-1))
    sum = 0.0
    nu = nuMin

    theta = boltzmann*tempe/(electronMass*cLight2)
    ne = eDens(r, tempi, tempe, v)
    h = height(r, tempi, tempe)
    temp_SSD = np.where(r > innerRadiusSSD*schwRadius, SSD_temp(r), 0.0)
    temp_SSD_inner = SSD_temp(innerRadiusSSD*2.0*schwRadius)
    k_es = ne*thomson
    tau_es = 2.0*k_es*h
    P = 1.0 - np.exp(-tau_es)
    A = 1.0 + 4.0*theta*(1.0+4.0*theta)
    kap = P*(A-1.0)/(1.0-P*A)
    phi = -np.log(P)/np.log(A)

    for i in range(iter):
        x = planck*nu/(electronMass*cLight2)/(3.0*theta)
        bnuSSD = bb_Intensity(nu, temp_SSD)
        bnuSSD_inner = bb_Intensity(nu, temp_SSD_inner) * np.square(r/(2.0*innerRadiusSSD*schwRadius))
        bnu = bb_Intensity(nu, tempe)
        bnu2 = 4.0*np.pi*bnu
        xi = xiBremss(nu, r, tempi, tempe, v) + xiSync(nu, r, tempi, tempe, v)
        k_nu = np.where(bnu2 > 1.0e-2*xi*h, xi/bnu2, 100.0)
        #l_eff = 1.0/np.sqrt(k_nu*(k_nu+k_es))
        #tau_es = 2.0*k_es*np.where(l_eff<h,l_eff,h)
        #P = 1.0-np.exp(-tau_es)
        #kap = P*(A-1.0)/(1.0-P*A)
        #phi = -np.log(P)/np.log(A)
        etaa = 1.0 + kap*(1.0-np.power(x, phi-1.0))
        etaaCD = np.where(etaa > 1.0, etaa - 1.0, 0.0)
        tau_nu = np.sqrt(np.pi)/2.0 * k_nu * h
        fluxxCD = np.where(SSDdisk == 1, np.where(innerRadiusSSD < 3.0, np.where(r < 2.0*innerRadiusSSD*schwRadius, bnuSSD_inner, bnuSSD), bnuSSD), 0.0)
        fluxx = np.where(bnu > 0.0 and tau_nu > 0.0, 2.0*np.pi/np.sqrt(3.0) * bnu * \
            (1.0-np.exp(-2.0*np.sqrt(3.0)*tau_nu)), 0.0)
        dnu = nu*(np.sqrt(paso)-1.0/np.sqrt(paso))
        sum += (fluxx*etaa + fluxxCD*etaaCD)*dnu
        #sum += fluxx*etaa*dnu
        nu = nu*paso
    return sum / (2*h)

def qminus_esin_fast(r, tempi, tempe, v):
    nuMin = 1.0e6
    nuMax = 1.0e21
    iter = 30
    paso = np.power(nuMax/nuMin, 1.0/float(iter-1))
    sum = 0.0
    nu = nuMin

    bb = BlackBody(tempe * u.K)
    theta = boltzmann*tempe/(electronMass*cLight2)
    ne = eDens(r, tempi, tempe, v)
    h = height(r, tempi, tempe)
    k_es = ne*thomson
    tau_es = 2.0*k_es*h
    A = 1.0 + 4.0*theta*(1.0+4.0*theta)
    for i in range(iter):
        x = planck*nu/(electronMass*cLight2)/(3.0*theta)
        bnu2 = 4.0*np.pi*bb(nu*u.Hz).value

        xi = xiBremss(nu, r, tempi, tempe, v)+xiSync(nu, r, tempi, tempe, v)
        #k_nu = np.where(bnu2 > 1.0e-2*xi*h,xi/bnu2,50.0)
        k_nu = xi/bnu2
        #l_eff = 1.0/np.sqrt(k_nu*(k_nu+k_es))
        #tau_es = 2.0*k_es*np.where(l_eff<h,l_eff,h)

        s_exp = tau_es*(tau_es+1.0)
        etaMax = 3.0*boltzmann*tempe/(planck*nu)
        jm = np.log(etaMax)/np.log(A)

        gamma1 = sp.gammainc(jm+1.0, A*s_exp)
        gamma2 = sp.gammainc(jm+1.0, s_exp)
        aux2 = s_exp*(A-1.0)

        #etaa = np.where(aux2<200.0,np.exp(aux2)*(1.0-gamma1)+etaMax*gamma2,etaMax*gamma2)
        etaa = np.exp(aux2)*(1.0-gamma1)+etaMax*gamma2

        tau_nu = np.sqrt(np.pi)/2.0 * k_nu * h
        fluxx = 2.0*np.pi/np.sqrt(3.0) * bnu * \
            (1.0-np.exp(-2.0*np.sqrt(3.0)*tau_nu))
        dnu = nu*(np.sqrt(paso)-1.0/np.sqrt(paso))
        sum += fluxx*etaa*dnu
        nu = nu*paso
    return sum*2.0/(np.sqrt(np.pi)*h)


def Den_r(logr, logTi, logTe, logv, j):

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    Te = np.exp(logTe)
    TeR = Te*eMMW
    te = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)

    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    machNumber2 = vnorm*vnorm / cs2norm
    invM2 = alpha*alpha/machNumber2
    invM2e = delta*invM2
    invM2i = (1.0-delta)*invM2
    Qvcs = machNumber2-1.0

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    aaux = a_aux(1.0/te) + te*aux_der(te)

    A = etai/2.0
    B = etae/2.0
    C = Qvcs
    E = etai*(0.5*beta*(3.0+etai)-invM2i)
    F = etae*(beta*0.5*etai-invM2i)
    G = beta*etai+invM2i
    I = etai*(0.5*beta*etae-invM2e)
    J = etae*(beta*(aaux+0.5*etae)-invM2e)
    K = beta*etae+invM2e

    return A*F*K - A*G*J + C*E*J - C*F*I + B*G*I - B*E*K


def Num1_r(logr, logTi, logTe, logv, j):

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    Te = np.exp(logTe)
    TeR = Te*eMMW
    te = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)

    rnorm = np.exp(logr)
    r = schwRadius*rnorm

    rho = massDensity(r, TiR, TeR, v)
    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    l = -alpha * cs2norm/vnorm * rnorm + j
    lK = keplAngMom(rnorm)
    machNumber2 = vnorm*vnorm / cs2norm
    invM2 = alpha*alpha/machNumber2
    invM2e = delta*invM2
    invM2i = (1.0-delta)*invM2
    Qvcs = machNumber2-1.0
    f1 = 1.5-s + 1.0/(1-1.0/rnorm)
    f2 = (l*l-lK*lK) / (rnorm*rnorm*cs2norm)
    f = f1+f2

    Pgas = beta*rho*cs2

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    Factor1 = r / (Pgas*v)
    Qie = qie_approx(r, TiR, TeR, v) * Factor1
    Qmin = qminus_esin(r, TiR, TeR, v) * Factor1

    Factor2 = 2.0*alpha*(j/(vnorm*rnorm))
    aaux = a_aux(1.0/te) + te*aux_der(te)

    B = etae/2.0
    C = Qvcs
    D = f
    F = etae*(beta*0.5*etai-invM2i)
    G = beta*etai+invM2i
    H = Factor2*(1.0-delta)-invM2i-beta*(etai*f1+Qie)
    J = etae*(beta*(aaux+0.5*etae)-invM2e)
    K = beta*etae+invM2e
    L = Factor2*delta-invM2e - beta*(etae*f1-Qie+Qmin)

    return B*G*L - B*H*K + F*D*K - F*C*L + J*C*H - J*D*G


def rhs_r(logr, y):

    logTi, logTe, logv, j = y

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    ti = boltzmann*TiR/(protonMass*cLight2)
    Te = np.exp(logTe)
    TeR = Te*eMMW
    te = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)
    rnorm = np.exp(logr)
    r = schwRadius*rnorm

    rho = massDensity(r, TiR, TeR, v)
    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    l = -alpha * cs2norm/vnorm * rnorm + j
    lK = np.power(rnorm, 1.5)/(rnorm-1.0)/np.sqrt(2.0)

    machNumber2 = vnorm*vnorm / cs2norm
    invM2 = alpha*alpha/machNumber2
    invM2e = delta*invM2
    invM2i = (1.0-delta)*invM2
    Qvcs = machNumber2-1.0
    f1 = 1.5-s + 1.0/(1-1.0/rnorm)
    f2 = (l*l-lK*lK) / (rnorm*rnorm*cs2norm)
    f = f1+f2

    Pgas = beta*rho*cs2

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    Factor1 = r / (Pgas*v)
    Qie = qie_approx(r, TiR, TeR, v) * Factor1
    Qmin = qminus_esin(r, TiR, TeR, v) * Factor1

    Factor2 = 2.0*alpha*(j/(vnorm*rnorm))
    ae = 3.0-6.0/(4.0+5.0*te) + te*30.0/np.square(4.0+5.0*te)
    ai = 3.0-6.0/(4.0+5.0*ti) + ti * 30.0/np.square(4.0+5.0*ti)

    A = etai/2.0
    B = etae/2.0
    C = Qvcs
    D = f
    E = etai*(beta*(ai+0.5*etai)-invM2i)
    F = etae*(beta*0.5*etai-invM2i)
    G = beta*etai+invM2i
    H = Factor2*(1.0-delta)-invM2i-beta*(etai*f1+Qie)
    I = etai*(beta*0.5*etae-invM2e)
    J = etae*(beta*(ae+0.5*etae)-invM2e)
    K = beta*etae+invM2e
    L = Factor2*delta-invM2e - beta*(etae*f1-Qie+Qmin)

    DEN = A*F*K - A*G*J + C*E*J - C*F*I + B*G*I - B*E*K
    N1 = B*G*L - B*H*K + F*D*K - F*C*L + J*C*H - J*D*G
    N2 = D*G*I - D*E*K + C*E*L - C*H*I + A*H*K - A*G*L
    N3 = D*E*J - D*F*I + B*H*I - B*E*L + A*F*L - A*H*J

    rhs_1 = N1/DEN
    rhs_2 = N2/DEN
    rhs_3 = N3/DEN
    rhs_4 = 0.0

    print(np.log10(r/schwRadius), N1, DEN)
    return np.asarray((rhs_1, rhs_2, rhs_3, rhs_4), dtype=float)


def Num_beta(logr, logTi, logTe, logv, j):

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    theta_i = boltzmann*TiR/(protonMass*cLight2)
    Te = np.exp(logTe)
    TeR = Te*eMMW
    theta_e = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)
    rnorm = np.exp(logr)
    r = schwRadius*rnorm

    rho = massDensity(r, TiR, TeR, v)
    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    l = -alpha * cs2norm/vnorm * rnorm + j
    lK = np.power(rnorm, 1.5)/(rnorm-1.0)/np.sqrt(2.0)

    machNumber2 = vnorm*vnorm / cs2norm
    bM2 = beta * machNumber2
    Qvcs = machNumber2-1.0
    alpha2 = alpha*alpha

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    term1 = dlogOmK_dlogr(rnorm) + s - 1.0
    term2 = 2.0*alpha*j/(rnorm*(-vnorm))

    tacc = r/(-v)
    pressure = rho*cs2
    factor = tacc / pressure
    Qie = qie_approx(r, TiR, TeR, v) * factor
    Qmin = qminus_dermer_fast(r, TiR, TeR, v) * factor
    #Qmin = qminus_esin_fast(r,TiR,TeR,v) * factor

    ae = 3.0-6.0/(4.0+5.0*theta_e) + theta_e * 30.0/np.square(4.0+5.0*theta_e)
    ai = 3.0-6.0/(4.0+5.0*theta_i) + theta_i * 30.0/np.square(4.0+5.0*theta_i)

    aux1 = alpha2 * (1.0-delta)

    B = etae * (bM2 * 0.5 * etai - aux1)
    C = etai * bM2 + aux1
    D = machNumber2 * (beta*etai * term1 - term2 * (1.0-delta) + Qie) - aux1

    E = etai * (bM2 * (ai+0.5) - alpha2)
    F = etae * (bM2 * (ae+0.5) - alpha2)
    G = bM2 + alpha2
    H = machNumber2 * (beta*term1 - term2 + Qmin) - alpha2

    J = 0.5*etae
    K = Qvcs
    L = -term1 + (l*l-lK*lK) / (rnorm*rnorm*cs2norm)

    return B*(G*L - H*K) + F*(D*K - C*L) + J*(C*H - D*G)


def Den_beta(logr, logTi, logTe, logv, j):

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    theta_i = boltzmann*TiR/(protonMass*cLight2)
    Te = np.exp(logTe)
    TeR = Te*eMMW
    theta_e = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)
    rnorm = np.exp(logr)
    r = schwRadius*rnorm

    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    l = -alpha * cs2norm/vnorm * rnorm + j
    lK = np.power(rnorm, 1.5)/(rnorm-1.0)/np.sqrt(2.0)

    machNumber2 = vnorm*vnorm / cs2norm
    bM2 = beta * machNumber2
    Qvcs = machNumber2-1.0
    alpha2 = alpha*alpha

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    ae = 3.0-6.0/(4.0+5.0*theta_e) + theta_e * 30.0/np.square(4.0+5.0*theta_e)
    ai = 3.0-6.0/(4.0+5.0*theta_i) + theta_i * 30.0/np.square(4.0+5.0*theta_i)

    aux1 = alpha2 * (1.0-delta)

    A = etai * (bM2 * (ai+0.5*etai) - aux1)
    B = etae * (bM2 * 0.5 * etai - aux1)
    C = etai * bM2 + aux1

    E = etai * (bM2 * (ai+0.5) - alpha2)
    F = etae * (bM2 * (ae+0.5) - alpha2)
    G = bM2 + alpha2

    I = 0.5*etai
    J = 0.5*etae
    K = Qvcs

    return A*(F*K - G*J) + C*(E*J - F*I) + B*(G*I - E*K)


def rhs_beta(logr, y):

    logTi, logTe, logv, j = y

    Ti = np.exp(logTi)
    TiR = Ti*iMMW
    theta_i = boltzmann*TiR/(protonMass*cLight2)
    Te = np.exp(logTe)
    TeR = Te*eMMW
    theta_e = boltzmann*TeR/(electronMass*cLight2)
    v = -np.exp(logv)
    rnorm = np.exp(logr)
    r = schwRadius*rnorm

    rho = massDensity(r, TiR, TeR, v)
    cs2 = sqrdSoundVel(TiR, TeR)

    vnorm = v/cLight
    cs2norm = cs2/cLight2

    l = -alpha * cs2norm/vnorm * rnorm + j
    lK = np.power(rnorm, 1.5)/(rnorm-1.0)/np.sqrt(2.0)

    machNumber2 = vnorm*vnorm / cs2norm
    bM2 = beta * machNumber2
    Qvcs = machNumber2-1.0
    alpha2 = alpha*alpha

    etai = Ti/(Ti+Te)
    etae = 1.0-etai

    term1 = dlogOmK_dlogr(rnorm) + s - 1.0
    term2 = 2.0*alpha*j/(rnorm*(-vnorm))

    tacc = r/(-v)
    pressure = rho*cs2
    factor = tacc / pressure
    Qie = qie_approx(r, TiR, TeR, v) * factor
    if TeR > 1.0e8:
        Qmin = qminus_dermer_fast(r, TiR, TeR, v) * factor
    else:
        Qmin = 0.0

    ae = 3.0-6.0/(4.0+5.0*theta_e) + theta_e * 30.0/np.square(4.0+5.0*theta_e)
    ai = 3.0-6.0/(4.0+5.0*theta_i) + theta_i * 30.0/np.square(4.0+5.0*theta_i)

    aux1 = alpha2 * (1.0-delta)

    A = etai * (bM2 * (ai+0.5*etai) - aux1)
    B = etae * (bM2 * 0.5 * etai - aux1)
    C = etai * bM2 + aux1
    D = machNumber2 * (beta*etai * term1 - term2 * (1.0-delta) + Qie) - aux1

    E = etai * (bM2 * (ai+0.5) - alpha2)
    F = etae * (bM2 * (ae+0.5) - alpha2)
    G = bM2 + alpha2
    H = machNumber2 * (beta*term1 - term2 + Qmin) - alpha2

    I = etai
    J = etae
    K = Qvcs
    L = -term1 + (l*l-lK*lK) / (rnorm*rnorm*cs2norm)

    DEN = A*(F*K - G*J) + C*(E*J - F*I) + B*(G*I - E*K)
    N1 = B*(G*L - H*K) + F*(D*K - C*L) + J*(C*H - D*G)
    N2 = D*(G*I - E*K) + C*(E*L - H*I) + A*(H*K - G*L)
    N3 = D*(E*J - F*I) + B*(H*I - E*L) + A*(F*L - H*J)

    rhs_1 = N1/DEN
    rhs_2 = N2/DEN
    rhs_3 = N3/DEN
    rhs_4 = 0.0

    print(np.log10(r/schwRadius), N1, DEN)
    return np.asarray((rhs_1, rhs_2, rhs_3, rhs_4), dtype=float)


def event(logr, y):
    logTi, logTe, logv, j = y
    return Den_r(logr, logTi, logTe, logv, j)


event.terminal = True


def bounds(log10j, logTiOut, logTeOut, logvOut):

    l_in = np.power(10.0, log10j)
    print("l_in = ", l_in)

    y0 = np.array([logTiOut, logTeOut, logvOut, l_in])
    solution = solve_ivp(rhs_r, (np.log(rOut), 0.5), y0,
                         method='LSODA', vectorized=False, events=event)

    print(solution.message)
    result = Num1_r(solution.t[-1], solution.y[0][-1],
                    solution.y[1][-1], solution.y[2][-1], l_in)
    print("RESULT = ", result)
    print("rSonic = ", np.exp(solution.t[-1]))
    print()

    return result


def event_beta(logr, y):
    logTi, logTe, logv, j = y
    return Den_beta(logr, logTi, logTe, logv, j)


event_beta.terminal = True


def bounds_beta(log10j, logTiOut, logTeOut, logvOut):

    l_in = np.power(10.0, log10j)
    print("l_in = ", l_in)

    y0 = np.array([logTiOut, logTeOut, logvOut, l_in])
    solution = solve_ivp(rhs_beta, (np.log(rOut), 0.5), y0,
                         method='LSODA', vectorized=False, events=event_beta)

    print(solution.message)
    result = Num_beta(solution.t[-1], solution.y[0]
                      [-1], solution.y[1][-1], solution.y[2][-1], l_in)
    print("RESULT = ", result)
    print("rSonic = ", np.exp(solution.t[-1]))
    print()

    return result
