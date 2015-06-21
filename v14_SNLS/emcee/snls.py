import pymc
import numpy as np
from matplotlib import pyplot
from scipy.stats import norm
from scipy import integrate
import emcee

# Some constants
numDims = 3;
problemSpaceBounds = np.array( [[60, 80], [0, 1], [0, 1]] );
lowestLogLiklVal = -1000;
# For experiments
DEBUG_MODE = False 

# Log Likelihood function for SN
################################################################################
def snlsLogLikl(evalPts, obsData, numObsToUse = None):

  # Define some constants
  LIGHT_SPEED = 299792.48;
  RESOLUTION = 100; # For 1D integration

  # Prelims
  numObs = obsData.shape[0];
  numPts = evalPts.shape[0];
  if numObsToUse is None:
    numObsToUse = numObs;

  # Decompose Data
  useObsData = obsData[0:numObsToUse, :];  
  z = useObsData[:, [0] ];
  obs = useObsData[:, [1] ];
  obsErr = useObsData[:, [2] ];

  # Create arrays for storing outputs
  logJointProbs = np.zeros( (numPts, 1) );
  lumMeans = np.zeros( (numPts, numObsToUse) );

  # Now iterate through each of the evalPts
  for i in range(numPts):
    H = evalPts[i, 0]
    OmegaM = evalPts[i, 1]
    OmegaL = evalPts[i, 2]
    if (H <= problemSpaceBounds[0, 0] or H >= problemSpaceBounds[0, 1] or
       OmegaM <= problemSpaceBounds[1, 0] or OmegaM >= problemSpaceBounds[1, 1] or
       OmegaL <= problemSpaceBounds[2, 0] or OmegaL >= problemSpaceBounds[2, 1]) :
      currLumMeans = -np.inf * np.ones( (numObsToUse, 1) ); obsLogLikl = -np.inf;
    else:
       
      currLumMeans = np.zeros( (numObsToUse, 1) );

      f = lambda t: 1/ np.sqrt( OmegaM * (1+t)**3 + OmegaL );
      for obsIter in range(numObsToUse):
        dLC = LIGHT_SPEED * (1 + z[obsIter, 0]) / H ;
        dLI = integrate.romberg(f, 0, z[obsIter, 0] );
        dL = dLC * dLI;
        currLumMeans[obsIter, 0] = 5 * np.log10(dL) + 25;
      obsLikl = norm.pdf(obs, currLumMeans, obsErr);
      obsLogLikl = sum( np.log(obsLikl) );
      if obsLogLikl < lowestLogLiklVal:
        obsLogLikl = lowestLogLiklVal;

    obsLogLikl = obsLogLikl;
    logJointProbs[i,0] = obsLogLikl;
    logJointProbs[i,0] = logJointProbs[i,0] * numObs / numObsToUse;
    lumMeans[[i], :] = currLumMeans.transpose();

  return logJointProbs, lumMeans


# Wrapper for normalized coordinates
################################################################################
def snlsNormLogLikl(normEvalPts, obsData):
  
  # First compute Normalized version
#   print normEvalPts
  unNormH = problemSpaceBounds[0,0] + \
    normEvalPts[:, [0]] * (problemSpaceBounds[0,1] - problemSpaceBounds[0,0] );
  unNormPts = normEvalPts;
  unNormPts[:,[0]] = unNormH;

  # Now compute the log likelihood
  return snlsLogLikl(unNormPts, obsData);


def lnProb(value, obsData):
  valueArr = np.array([value]);
  ljps, lumMeans = snlsNormLogLikl(valueArr, obsData);
  ljps = float(ljps[0,0]);
  return ljps


# Unit test for evalSNLSLogLikl
################################################################################
def unittest_evalSNLSLogLikl(data):
  pts = np.array([ [68.0, 0.24, 0.68],
                   [61, 0.5, 0.5],
                   [65, 0.3, 0.65] ]);
  logJointProbs, lumMeans = snlsLogLikl(pts, data);
  print pts
  print logJointProbs.shape, lumMeans.shape, '\n', logJointProbs, 

  print "Test on Noramlized coords:"
  pts = np.array([ [0.4, 0.24, 0.68],
                   [0, 0.5, 0.5],
                   [0, 0, 0],
                   [0.1, 0.1, 0.1],
                   [0.5, 0.5, -1],
                   [0.25, 0.3, 0.65] ]);
  logJointProbs, lumMeans = snlsNormLogLikl(pts, data);
  print logJointProbs.shape, lumMeans.shape, '\n', logJointProbs, 


# Pymc set up for this likelihood function
################################################################################
def getPymcVariable(data):

  def logp(value):
    valarr = np.array([value]);
    ljps, lumMeans = snlsNormLogLikl(valarr, data);
    ljps = float(ljps[0,0]);
#     print type(ljps), str(value) + ' ' + str(ljps)
    return ljps

  def random():
    return np.random.rand(numDims);

  dtype = type(random());
  initPt = [0.45, 0.24, 0.68];
#   initPt = [0.5, 0.5, 0.5];

  ret = pymc.Stochastic(logp = logp,
                        doc = 'SNLS RV',
                        name = 'SNLS',
                        parents = {},
                        random = random,
                        trace = True,
                        value = initPt,
                        dtype = dtype,
                        observed = False,
                        cache_depth = 2,
                        plot = True,
                        verbose = 0 );
  return ret


# Main loop
################################################################################
if __name__ == '__main__':

  data = np.loadtxt('../davisdata.txt');
  unittest_evalSNLSLogLikl(data);

  numExperiments = 30;
  numWalkers = 6;
  numStops = 43;
  numPtsPerStop = 25;

  logLiklHandle = lambda x : lnProb(x, data);

#   sampler = emcee.EnsembleSampler(numWalkers, numDims, logLiklHandle);
#   numIters = np.ceil(numSamples/numWalkers);
#   print "Running for " + str(numIters) + " iterations" 
#   p0 = [np.random.rand(numDims) for i in xrange(numWalkers)];
#   sampler.run_mcmc(p0, numIters);
#   print "Mean Acceptance fRaction: " + str(np.mean(sampler.acceptance_fraction));
#   allSamples = sampler.flatchain[:,0];
#   print allSamples.shape
#   pyplot.hist(allSamples, 100)
#   pyplot.show()

  for expIter in range(numExperiments):

    # Random starting 
    print 'Experiment ' + str(expIter)
    p0 = [np.random.rand(numDims) for i in xrange(numWalkers)];
    sampler = emcee.EnsembleSampler(numWalkers, numDims, logLiklHandle);
    pos = p0;

    for stopIter in range(numStops):

      fileName = 'emceeRes/exp' + str(expIter) + '_stop' + str(stopIter) + '.txt'
      pos, prob, state = sampler.run_mcmc(p0, numPtsPerStop);
      allSamples = sampler.flatchain;
      print '  Stop ' + str(stopIter) + ', ' + str(allSamples.shape)
      np.savetxt(fileName, allSamples); 

#     fileName = 'experiment' + str(expIter) + '.txt';
#     SNLS = getPymcVariable(data)
#     model = pymc.Model([SNLS]);
#     mcmcModel = pymc.MCMC(model);
# #     mcmcModel.use_step_method(pymc.Metropolis, SNLS, proposal_sd=1);
#     mcmcModel.sample(numSamples, burn=0, thin=1);
#     x = mcmcModel.trace('SNLS')[:]
#     np.savetxt(fileName, x);
    
#     H = x[:,0];
#     OM = x[:,1];
#     OL = x[:,2];
#     pyplot.figure()
#     pyplot.hist( H, bins=30);
#     pyplot.figure()
#     pyplot.hist( OM, bins=30);
#     pyplot.figure()
#     pyplot.hist( OL, bins=30);
#     pyplot.show();
    
#     print x.shape, type(x),

