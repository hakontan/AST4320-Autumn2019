import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import astropy.units as unit




T_CMB = 2.725*unit.K
a = np.linspace(1e-4, 1, 1000)
z = 1/a - 1

def T_rad(a):
    return(T_CMB * np.power(a, -1))

def T_gas(a):
    return(T_CMB/1091 * np.power(a, -2))

plt.loglog(a, T_rad(a), label="T_rad")
plt.loglog(a, T_gas(a), label="T_gas")
plt.xlabel("log(a)"); plt.ylabel("log(T)")
plt.legend()
plt.savefig("exercise3.eps")
plt.show()