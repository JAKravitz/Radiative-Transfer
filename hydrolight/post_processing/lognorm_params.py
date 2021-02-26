# for randomly selecting out of lognormal dist.

def lognorm_params(mode, stddev):
    import numpy as np
    """
    Given the mode and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    p = np.poly1d([1, -1, 0, 0, -(stddev/mode)**2])
    r = p.roots
    sol = r[(r.imag == 0) & (r.real > 0)].real
    shape = np.sqrt(np.log(sol))
    scale = mode * sol
    return shape, scale

def lognorm_random(sigma,scale,num):
    import numpy as np
    mu = np.log(scale)
    s = np.random.lognormal(mu, sigma, num)
    return s