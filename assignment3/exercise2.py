import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import astropy.constants as const


c        = const.c.cgs            # speed of light in cgs
H0       = cosmo.H(0).cgs         # Hubble constant in cgs
Sigma_T  = 6.65e-25 * u.cm**2     # Thomson cross section of molecule [cm^2]
O_m      = 0.308                  # energy density parameter for mass 
O_lambda = 0.692                  # energy density parameter for lambda


def tau_e(z):
    """
    Returns the optical depth for the ionized intergalactic medium
    as a function of redshift z.

    Parameters:
    -----------
    z: redshift at which the optical depth is to be calculated
    """
    _z = np.linspace(0, z, 100) # New redshift array to define limits in the integral
    numerator = Sigma_T * 1.9e-7 * u.cm**-3 * (1 + _z )**2
    denominator = H0 * np.sqrt(O_m * (1 + _z)**3 + O_lambda)
    integral = c * numerator / denominator
    res = np.trapz(integral, _z)
    return(res)


z = np.linspace(0, 10, 100) # Redshift array
tau = np.zeros(len(z))
for i in range(len(z)):
    tau[i] = tau_e(z[i])

fig, ax = plt.subplots()
ax.plot(z, tau)
ax.set_xlabel("Redshift (z)")
ax.set_ylabel(r"$\tau_{e}(z)$")
fig.savefig("exercise2.pdf")
plt.show()