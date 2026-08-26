import numpy as np
import scipy.constants as const
import json
import os

# CONSTANTS [CGS]

G = const.gravitational_constant * 1.0e3
cLight = const.speed_of_light * 1.0e2
cLight2 = cLight*cLight
solarMass = 1.98847e33
eddAccRate = 1.39e18

def virialTemp(r):
    return 3.6e12/r

adafFile = 'adafFile.txt'
adafParameters = 'adafParameters.txt'

# Fix parameters
coronaFactor = 1.0
SSDdisk = 0
pIndex = 0.1
correctorAccRate = 0.01

# INITIAL PARAMETERS

blackHoleMass = 6.5e9       # [Msol]
eddAccRate = eddAccRate * blackHoleMass
accRateNorm = 0.00042          # [MdotEdd]
rOut = 1e4                  # [Rs]
accRateOut = accRateNorm * eddAccRate * coronaFactor
accRateCD = accRateNorm * eddAccRate
innerRadiusSSD = 0.9e4        # [Rs]

# Adimensional Disk parameters
beta = 0.9
alpha = 0.1
delta = 0.3
s = 0.35

# EIGENVALUES
log10j0 = np.log10(1.2163)
log10j1 = np.log10(2.1)

# Boundary conditions

# rOut ~ 10^2 rS
#temp_i_Out = 0.6*virialTemp(rOut)
#temp_e_Out = 0.08*virialTemp(rOut)
#lamda = 0.5

if np.abs(np.log10(rOut/1.0e2)) < np.abs(np.log10(rOut/1.0e4)):
    temp_i_Out = 0.6*virialTemp(rOut)
    temp_e_Out = 0.08*virialTemp(rOut)
    lamda = 0.5
else:
    temp_i_Out = 0.2*virialTemp(rOut)
    temp_e_Out = 0.19*virialTemp(rOut)
    lamda = 0.2

# The historical scripts remain directly runnable, but the automated pipeline can
# supply the same inputs without editing this file.  Values in RIAF_CONFIG_JSON
# override the defaults above before global_eqs imports them.
_config_path = os.environ.get("RIAF_CONFIG_JSON")
if _config_path:
    with open(_config_path, encoding="utf-8") as _config_file:
        _cfg = json.load(_config_file)
    for _name in (
        "blackHoleMass", "accRateNorm", "rOut", "innerRadiusSSD",
        "beta", "alpha", "delta", "s", "pIndex", "coronaFactor",
        "SSDdisk", "correctorAccRate", "log10j0", "log10j1",
        "temp_i_Out", "temp_e_Out", "lamda",
    ):
        if _name in _cfg:
            globals()[_name] = _cfg[_name]
    eddAccRate = 1.39e18 * blackHoleMass
    accRateOut = accRateNorm * eddAccRate * coronaFactor
    accRateCD = accRateNorm * eddAccRate
    if "temp_i_Out" not in _cfg or "temp_e_Out" not in _cfg or "lamda" not in _cfg:
        if np.abs(np.log10(rOut/1.0e2)) < np.abs(np.log10(rOut/1.0e4)):
            temp_i_Out = 0.6*virialTemp(rOut)
            temp_e_Out = 0.08*virialTemp(rOut)
            lamda = 0.5
        else:
            temp_i_Out = 0.2*virialTemp(rOut)
            temp_e_Out = 0.19*virialTemp(rOut)
            lamda = 0.2
