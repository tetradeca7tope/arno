import cern.jet.stat.Probability;
import ParamVector;


/** Computes the chisq distance and alpha value for Tegmark's LSS data. */
public class LssModeller {
  private final double dof = 14.0; // 20 parameters, 6 dof.
  LssWrapper wrap = new LssWrapper();
  
  /** Returns chisq distance */
  public double computeValue(ParamVector vector) {
    return computeValue(vector.getLssVector());
  }
  
  public double computeValue(double[] parameters) {
    return wrap.computeChisq(parameters);
  }

  public double getAlpha(double dist) {
    return 1.0-Probability.gamma(0.5, dof/2.0, dist);
  }
}

