import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
from astropy.cosmology import WMAP9 as cosmo

N = 10000
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

def solver(Om, Ol, plot=False):
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
    if plot==True:
        plt.loglog(a, delta, label = r"$\Omega _m$ = {0}, $\Omega _m$ = {1}"
                                        .format(Om, Ol))
        plt.xlabel("log(a)"); plt.ylabel(r"log($\delta$)")

    return(delta, deltadot)

def growth_factor(Om, Ol):
        z = 1/a - 1
        delta, deltadot = solver(Om, Ol)
        f = np.power(Hubble(a, Om, Ol), -1) * (deltadot/delta)
        plt.loglog(z, f, label = r"$\Omega _m$ = {0}, $\Omega _m$ = {1}"
                                    .format(Om, Ol))
        plt.xlabel("log(z)")
        plt.ylabel(r"log($F=\frac{dln\delta}{dlna}$)")
        


solver(1, 0, plot=True)
solver(0.3, 0.7, plot=True)
solver(0.8, 0.2, plot=True)
plt.legend()
plt.title("Overdensity")
plt.savefig("exercise2dot2.eps")
plt.show()

growth_factor(1, 0)
growth_factor(0.3, 0.7)
growth_factor(0.8, 0.2)
plt.legend()
plt.title("Growth factor")
plt.savefig("exercise2dot3.eps")
plt.show()