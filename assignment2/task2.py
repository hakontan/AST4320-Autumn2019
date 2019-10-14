import numpy as np
import matplotlib.pyplot as plt
from progress.bar import IncrementalBar 

N = int(1e3)


def variance(Sc):
    """
    Returns the variance as a function of smoothing scale.
    """
    return(np.pi / Sc**4)

def S(var):
    """
    calculating the smoothing coefficient.
    """
    return((np.pi / var)**(1/4))


def P(S, delta):
    """
    Probability distribution for regular
    random walk.
    """
    var = variance(S)
    term1 = 1 / np.sqrt(2*np.pi*var) 
    term2 = np.exp(-delta**2 / (2*var))
    return(term1 * term2)

def P_noncross(S, delta):
    """
    Probability distribution for random walk
    where delta does not cross delta_crit.
    """
    var = variance(S)
    term1 = 1 / np.sqrt(2 * np.pi * var)
    term2 = np.exp(- 0.5 * delta**2 / var )
    term3 = np.exp(- 0.5 * (2 - delta)**2 / var)
    return(term1 * (term2 - term3))


def randomwalk(N):
    """
    Computes the regular random walk.
    """
    delta_arr = np.zeros(N)

    bar = IncrementalBar("Processing", max = N)
    for i in range(N):
        eps = 5e-3
        
    
        var_1 = 1e-4
        Sc = S(var_1)
        delta = np.random.normal(loc=0, scale=np.sqrt(var_1), size=1)

        while Sc >= 1:
            Sc -= eps
            var_2 = variance(Sc)
            var_12 = var_2 - var_1
            beta = np.random.normal(loc=0.0, scale=np.sqrt(var_12), size=1)
            delta += beta
            var_1 = var_2
        bar.next()
        delta_arr[i] = delta

    bar.finish()

    delta_sort = np.sort(delta_arr)
    P_dist = P(1, delta_sort)

    plt.plot(delta_sort, P_dist, color="b")
    plt.hist(delta_arr, bins="auto", density=True)
    plt.savefig("randomwalk.eps")
    plt.show()


def randomwalk_noncross(N):
    """
    Computes the random walk where delta does not cross
    delta_crit.
    """
    bar = IncrementalBar("Processing", max = N)
    delta_arr = np.zeros(N)
    for i in range(N):
        eps = 5e-3
        

        var_1 = 1e-4
        Sc = S(var_1)
        delta = np.random.normal(loc=0.0, scale=np.sqrt(var_1), size=1)

        while Sc >= 1:
            Sc -= eps
            var_2 = variance(Sc)
            var_12 = var_2 - var_1
            beta = np.random.normal(loc=0.0, scale=np.sqrt(var_12), size=1)
            if delta >= 1:
                delta_arr[i] = np.NaN
                break
            delta += beta
            var_1 = var_2

        if delta_arr[i] == np.NaN:
            continue
        elif delta_arr[i] == 0:
            delta_arr[i] = delta


        bar.next()

    bar.finish()

    delta_arr_real = delta_arr[np.where(np.isnan(delta_arr) == False)]
    delta_sort = np.sort(delta_arr_real)
    
    P_dist = P_noncross(1, delta_sort)
    #print(P_dist)
    integral = np.trapz(y = P_dist, x = delta_sort)
    
    plt.plot(delta_sort[np.where(P_dist >= 0)], P_dist[np.where(P_dist >= 0)]/integral, color="b")
    plt.hist(delta_sort, bins="auto", density=True)
    plt.savefig("randomwalk_noncross.eps")
    plt.show()

randomwalk(N)
randomwalk_noncross(N)