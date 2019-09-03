import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
from astropy.cosmology import WMAP9 as cosmo

N = 1000
H0 = cosmo.H0

a0 = 1e-3
a = np.linspace(a0, 1, N)
da = a[1] - a[0]

def D_deltadot(H, a , delta, deltadot):
    term1 = 3*H*delta / (2*a)
    term2 = 2*deltadot / a
    return(term1 - term2)

def D_delta(deltadot, adot):
    return(deltadot / adot)

def Hubble(a, Om, Ol):
    return(H0.value * np.sqrt(Om*a**(-3) + Ol))

def integrate(Om, Ol):
    delta = np.zeros(N)
    deltadot = np.zeros(N)

    delta[0] = 1e-3
    deltadot[0] = a0*Hubble(a0, Om, Om)
    
    for i in range(N-1):
        H = Hubble(a[i], Om, Ol)

        deltadot[i+1] = (deltadot[i] + da * D_deltadot(H,
                                                       a[i],
                                                       delta[i],
                                                       delta[i+1]))

        delta[i+1]    = (delta[i] + da * D_delta(deltadot[i+1],
                                                 a[i] * H))
    plt.loglog(a, delta, label = r"$\Omega _m$ = {0}, $\Omega _m$ = {1}"
                                   .format(Om, Ol))
    
integrate(1, 0)
integrate(0.3, 0.7)
integrate(0.8, 0.2)

plt.legend()
plt.show()