class LssWrapper {
  
  static {
    System.loadLibrary("lss_wrapper");
  }
  
  public LssWrapper() {
    System.out.println("Initializing Fortran.");
    initialize();
  }

  private native void initialize();
  
  private native double computeChisqC(double[] array);
  
  /** Computes the chisq values of the lss data against the given parameter
 * vector.
   * Note: this code uses an 6D parameter vector, not the typical 7D vector! */
  public double computeChisq(double[] params) {
    if (params.length != 9) {
      throw new IllegalArgumentException("Invalid Number of arguments");
    }

    /*
    double[] x = new double[9];

    x[0] = params[0];     // Omega_k 
    x[1] = params[1];     // Omega_DE
    x[2] = params[2];     // omega_c
    x[3] = params[3];     // omega_b
    x[4] = params[4];     // n_s
    x[5] = 0.6845;        // Ignore: A_s
    x[6] = 0.0;           // Ignore: alpha
    x[7] = params[5];     // galaxy bias 
    x[8] = 30.81;         // Ignore: Q_nl
    */
    
    return computeChisqC(params);
  }
  
  
  public static void main(String[] args) {
    System.out.println("Calling Init");
    LssWrapper w = new LssWrapper();

    // best fit model.
    double[] x = {0.0, 0.762, 0.1045, 0.02233, 0.951, 0.6845, 0.0, 1.908,
30.81};
    double chisq = w.computeChisq(x);
    System.out.println("Chisq: "+chisq);
  }
}

