import pymc
import numpy as np
from matplotlib import pyplot

numDims = 2;


# The pymc stochastic RV
def getGaussPYMCObject( mu, sigma, initPt = None ):

  numDims = len(mu);
  
  if initPt is None:
    initPt = mu;

  def logp(value):
    ret = 0;
    for d in range(numDims):
      ret += np.sum( -np.log(sigma[d]) - 0.5 * np.log(2 * np.pi) -
                   (value[d]-mu[d])**2 / 2*(sigma[d]**2) ) - 25;
    return ret

  def random():
    return np.random.rand(numDims);

  dtype = type(mu);
  print dtype,

  result = pymc.Stochastic( logp = logp,
                            doc = 'Gaussian RV',
                            name = 'X',
                            parents = {},
                            random = random,
                            trace = True,
                            value = initPt,
                            dtype = dtype,
                            observed = False,
                            cache_depth = 2,
                            plot = True,
                            verbose = 0 );
  return result


# Main
if __name__ == '__main__':
  
  mu = [0.5, 0.5];
  sigma = [1.4, 2.3];
  initPt = [0, 0];
#   mu = np.array([[1.0, 2.0]]);
#   sigma = np.array([[1.4, 2.3]]);
#   initPt = np.array([[-10.1, 10.3]]);
#   mu = np.array([[1.0]]);
#   sigma = np.array([[1.4]]);
#   initPt = np.array([[-10.1]]);
#   mu = np.array([[1.0], [2.0]]);
#   sigma = np.array([[1.4], [2.3]]);
#   initPt = np.array([[-10.1], [10.3]]);
  X = getGaussPYMCObject(mu, sigma, initPt);

  numSamples = 5000;
  model = pymc.Model([X]);
  mcmcModel = pymc.MCMC(model);
  mcmcModel.use_step_method(pymc.Metropolis, X, proposal_sd = 1);
  mcmcModel.sample(numSamples, burn=0, thin=1);
  x = mcmcModel.trace('X')[:]
  print x.shape

  x0 = x[:, 0].reshape(numSamples, 1);
  x1 = x[:, 1].reshape(numSamples, 1);
  pyplot.figure()
  pyplot.hist(x0, bins=30)
  pyplot.figure()
  pyplot.hist(x1, bins=30)
  pyplot.show();
  
