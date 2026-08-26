# Physical and astrophysical constants in CGS

import numpy as np
import astropy.constants as const
import scipy.constants as const_scpy

# Physical constants

cLight = const.c.cgs.value					# Speed of light in vacuum
cLight2 = cLight*cLight						# Speed of light in vacuum squared
electronMass = const.m_e.cgs.value			# Electron mass
electronRestEnergy = electronMass*cLight2	# Electron rest energy
protonMass = const.m_p.cgs.value			# Proton mass
protonRestEnergy = protonMass*cLight2		# Proton rest energy
planck = const.h.cgs.value					# Planck constant
planckbar = const.h.cgs.value/(2*np.pi)					# Planck bar constant (h partida)
boltzmann = const.k_B.cgs.value				# Boltzmann constant
electronCharge = const.e.gauss.value		# Electron charge
electronRadius =  const_scpy.physical_constants['classical electron radius'][0]*100		# Classical electron radius
thomson = const.sigma_T.cgs.value			# Thomson cross section
fineStructConst = 1.0/137.0					# Fine structure constant
gravConstant = const.G.cgs.value			# Gravitational constant
stefanBoltzmann = const.sigma_sb.cgs.value	# Stefan-Boltzmann constant
atomicMassUnit = const.u.cgs.value			# Atomic mass unit

# Astrophysical constants

solarMass = const.M_sun.cgs.value			# Solar mass
parsec = const.pc.cgs.value					# 1 parsec
eddLuminosity_1Msol = 1.26e38
eddAccRate_1Msol = 10 * eddLuminosity_1Msol / cLight2

# Unit conversion

Jy = 1.0e-23                    # 1 Jansky
eV = 1.602e-12					# 1 electronVolt
yrs_to_sec = 365.25 * 24 * 3600
Hz_to_eV = 1/2.417989242e14 

# Geometrized units conversion

