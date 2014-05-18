package edu.cmu.cs.auton.cosmo.lss;

import java.io.*;
import java.util.Arrays;

import edu.cmu.cs.auton.cosmo.AbstractPlotter;
import edu.cmu.cs.auton.spatial.objects.HyperRectangle;
import edu.cmu.cs.auton.util.RangeEntry;
import edu.cmu.cs.auton.util.RangeIterator;
import edu.cmu.cs.auton.util.StopWatch;

public class LssPlotter extends AbstractPlotter {
  final static LssModeller modeller = new LssModeller();
  
  private final static String[] labels =  {"{/Symbol=8 W}_k", "{/Symbol=8
W}_{DE}", 
    "{/Symbol=8 w}_{c}", "{/Symbol=8 w}_{b}", "Spectral Index", "Galaxy Bias"};  
    
  private final static HyperRectangle paramRanges = new HyperRectangle(new
double[]{
      0.0, 1.0,    // Omega_k
      0.0, 1.0,    // Omega_DE
      0.0, 3.0,    // omega_c
      0.001, 0.25, // omega_b
      0.5, 1.7,    // n_s
      0.1, 3.0     // galaxy bias
  });

  private double[][][] views;
  private StopWatch watch = new StopWatch();

  public LssPlotter() {
    super(paramRanges, labels);
  }
  
  
  public String getPrefix() {
    return "tegmarkLss";
  }
  
  private double[] getVector(double[] x) {
    double[] result = new double[9];
    result[0] = x[0];
    result[1] = x[1];
    result[2] = x[2];
    result[3] = x[3];
    result[4] = x[4];
    result[5] = 0.6845; // Ignore: A_s
    result[6] = 0.0;    // Ignore: alpha
    result[7] = x[5];
    result[8] = 30.81;  // Ignore: Q_nl
    return result;
  }
  
  public void createPlots(int numPtsPerDim) {
    System.out.println("NumPtsPerDim: "+numPtsPerDim);
    this.views = new double[combinations.length][numPtsPerDim][numPtsPerDim];

    RangeIterator iter = new RangeIterator(ranges, numPtsPerDim);

    double cutoff = Math.exp(12.59159/-2.0);
    int numLines = 0;
    
    try {
      PrintWriter goodOut = new PrintWriter(new FileWriter(new
File(getPrefix()+".dat")));
      
      watch.start();
      while (iter.hasNext()) {
        RangeEntry e = iter.next();
        int[] indexes = e.getIndexes();
        double[] x = e.getPoint();
        
        
        
        double distance = modeller.computeValue(getVector(x));
        double alpha = modeller.getAlpha(distance);
        
        // System.out.println(Arrays.toString(x)+" "+distance);
        
        if (!Double.isNaN(alpha) && alpha >= cutoff) {
          goodOut.printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f    %15.10f
%15.10f%n",
              x[0], x[1], x[2], x[3], x[4], x[5], distance, alpha);
          
          for (int i=0; i<combinations.length; i++) {
            int xIndex = indexes[combinations[i][0]];
            int yIndex = indexes[combinations[i][1]];
            
            if (views[i][xIndex][yIndex] < alpha) {
              views[i][xIndex][yIndex] = alpha;
            }
          }
        }
        
        if (e.isNewLine()) {
          numLines++;
          
          if ((numLines % 100) == 0) {
            System.out.println(numLines+"  "+Arrays.toString(indexes)+"
"+watch);
          }
        }
      }
      goodOut.close();
      watch.end("Iterations Complete: ");

      write(views);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  public static void main(String[] args) {
    LssPlotter plotter = new LssPlotter();
    plotter.createPlots(30);
  }
}

