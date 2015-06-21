import numpy as np
import emcee

def lnprob(x, ivar):
  return -0.5 * np.sum(ivar * x**2)

if __name__ == '__main__':
  ndim, nwalkers = 10, 2;
  ivar = 1./ np.random.rand(ndim);
  p0 = [np.random.rand(ndim) for i in range(nwalkers)];
  print ivar
  print p0

  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar]);
  sampler.run_mcmc(p0, 1000);
