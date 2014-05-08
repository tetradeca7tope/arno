import numpy as np
import emcee


def evalSNSLLogLikl(evalPts, obsData, numObsToUse):

  LIGHT_SPEED = 299792.458; # in Km/s
  RESOLUTION = 100; # resolution for 1D integration

  # Prelims
  (numPts, temp) = evalPts.size;
  (numObs, temp) = obsData.size;

  # Decompose obsData
  z = obsData


  return 0


if __name__ == '__main__':
  
  
