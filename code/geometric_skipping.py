from __future__ import division
import math
import numpy as np
import numpy.random
#from scipy.stats import geom
from cryptorandom.cryptorandom import SHA256


def geometric_skipping(population, p, seed):
    if isinstance(seed, int):
        ss = SHA256(seed)
    else:
        raise ValueError('%r cannot be used to seed SHA256 instance' % seed)
    
    n = len(population)
    sample = []
    m = int(n*p)
#    rvs = geom.rvs(p=p, size=m) # this is what we'd do if we were using numpy.random
    runifs = ss.random(m)
    rvs = list(map(lambda u: int(1 + np.log(u)/np.log(1-p)), runifs))
    val = np.cumsum(rvs)-1
    for i in range(m):
        if val[i] < n:
            sample.append(population[val[i]])
        else:
            return sample
    index = val[-1]
    while index < (n-1):
        index += int(1 + np.log(ss.random())/np.log(1-p))
        if index < n:
            sample.append(population[index])
    return sample
    


def test_transformation():
    u = np.array(range(100))
    runifs = u/len(u)
    np.testing.assert_equal(runifs[0], 0)
    p = 0.5
    rvs = list(map(lambda u: int(1 + np.log(u)/np.log(1-p)), runifs[1:]))
    np.testing.assert_equal(rvs[0], 7) # floor( 1 + ln(0.01)/ln(0.5) )
    np.testing.assert_equal(rvs[-1], 1) # floor( 1 + ln(0.99)/ln(0.5) )


def test_geometric_generation():
    ss = SHA256(12345)
    runifs = ss.random(10000)
    p = 0.5
    rvs = list(map(lambda u: int(1 + np.log(u)/np.log(1-p)), runifs))
    np.testing.assert_almost_equal(np.mean(rvs), 1/p, 1) # expected value of geometric(p)
    
    
def test_geometric_skipping():
    population = list(range(100))
    rvs = [5, 5, 10, 10, 20, 20, 30, 40]
    val = np.cumsum(rvs)-1 # indices to be chosen. array([  4,   9,  19,  29,  49,  69,  99, 139])
    n = len(population)
    m = len(val)
    sample = []
    for i in range(m):
        if val[i] < n:
            sample.append(population[val[i]])
        else:
            break
    np.testing.assert_equal(sample, val[:-1])
    
    sample = []
    rvs = [5, 5, 10, 10, 20, 20]
    p = 0.5
    val = np.cumsum(rvs)-1 # indices to be chosen. array([  4,   9,  19,  29,  49,  69])
    m = len(val)
    for i in range(m):
        if val[i] < n:
            sample.append(population[val[i]])
        else:
            break
    index = val[-1]
    while index < (n-1):
        index += int(1 + np.log(0.5)/np.log(1-p)) # add 2 to index
        if index < n:
            sample.append(population[index])
    np.testing.assert_equal(sample, list(val) + list(range(71, 100, 2)))

if __name__ == "__main__":
    test_transformation()
    test_geometric_generation()
    test_geometric_skipping()
