public class ParamVector {
  public enum Param {tau, OmegaDE, OmegaM, omegaDM, omegaB, fn, ns, 
    galaxy_bias, OmegaB, OmegaC, OmegaN, H0, OmegaK, omegaC, OmegaT}; 
    

  
  /** store 8D representation: 
   * tau, OmegaDE, OmegaM, omegaDM, omegaB, fn, ns, galaxy_bias */
  final double[] values;
  
  public ParamVector(double[] values) {
    this.values = values;
  }
  
  public double getParam(Param p) {
    return getParam(values, p);
  }

  public double[] getValues() {
    return values;
  }

  public double[] getSupernovaVector() {
    return getSupernovaVector(values);
  }
  
  public double[] getWmapVector() {
    return getWmapVector(values);
  }
  
  public double[] getLssVector() {
    return getLssVector(values);
  }
  
  
  // ------------------------------ Static Methods --------------------------
  public static double[] getSupernovaVector(double[] values) {
    // supernova vector is: H0, x.OmegaM, x.OmegaDE

    double[] x = new double[3];
    
    x[0] = getParam(values, Param.H0);
    x[1] = getParam(values, Param.OmegaM);
    x[2] = getParam(values, Param.OmegaDE);
    
    return x;
  }
  
  public static double[] getWmapVector(double[] values) {
    double[] x = new double[7];
    System.arraycopy(values, 0, x, 0, 7);

    return x;
  }
  
  public static double[] getLssVector(double[] values) {
    // lss vector is: OmegaK, OmegaDE, omegaC, omegaB, ns, As, alpha, galaxy
    // bias, Q_nl 
    
    double[] x = new double[9];
    x[0] = getParam(values, Param.OmegaK);
    x[1] = getParam(values, Param.OmegaDE);
    x[2] = getParam(values, Param.omegaC);
    x[3] = getParam(values, Param.omegaB);
    x[4] = getParam(values, Param.ns);
    x[5] = 0.6845; // Ignore: A_s
    x[6] = 0.0;    // Ignore: alpha
    x[7] = getParam(values, Param.galaxy_bias);
    x[8] = 30.81;  // Ignore: Q_nl

    return x;
  }
  
  public static double getParam(double[] values, Param p) {
    switch (p) {
    case tau: 
      return values[0];
      
    case OmegaDE:
      return values[1];
      
    case OmegaM: 
      return values[2];
      
    case omegaDM:
      return values[3];
      
    case omegaB:
      return values[4];
      
    case fn:
      return values[5];
      
    case ns:
      return values[6];
      
    case galaxy_bias:
      return values[7];
      
    case OmegaB: {
      double h = geth(values);
      return values[4]/(h*h);  // OmegaB
    }
    
    case OmegaC: {
      double h = geth(values);
      return values[3]/(h*h) *(1.0-values[5]); //OmegaC   
    }

    case omegaC:
      return values[3]*(1.0-values[5]); //omegaC   
    
    case OmegaN: {
      double h = geth(values);
      return values[3]/(h*h)*values[5];   // OmegaN
    }
    
    case H0:
      return 100.0*geth(values);

    case OmegaK:
      return 1.0 - values[2] - values[1];
      
    case OmegaT: 
      return values[1]+values[2];
      
    default: 
      throw new RuntimeException();
    }
  }
  
  private static double geth(double[] values) {
    // v3 = omegaDM, v4= omegaB, v2 = OmegaM 
    return Math.sqrt((values[3]+values[4])/values[2]);
  }
}

