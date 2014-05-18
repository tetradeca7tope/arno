package edu.cmu.cs.auton.cosmo.lss;

import edu.cmu.cs.auton.cosmo.AbstractIteratorPlotter;
import edu.cmu.cs.auton.cosmo.ParamVector;

public class LssPlotter2 extends AbstractIteratorPlotter {
  final static LssModeller modeller = new LssModeller();
  
  public LssPlotter2() {
    super(defaultRanges,
        defaultLabels,
        new int[]{1,2,3,4,5,6,7});
  }
  
  @Override
  protected double computeValue(double[] x) {
    double distance = modeller.computeValue(ParamVector.getLssVector(x));
    double alpha = modeller.getAlpha(distance);
    return alpha;
  }

  @Override
  public String getPrefix() {
    return "lss.new";
  }
  
  public static void main(String[] args) {
    LssPlotter2 plotter = new LssPlotter2();
    plotter.createPlot(20);
  }
}

