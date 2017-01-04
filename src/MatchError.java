public class MatchError {
  public static boolean almostEquals(double[] d1, double[] d2, double eps) {
    for (int i = 0; i < d1.length; i++) {
      double v1 = d1[i];
      double v2 = d2[i];
      if (!almostEquals(v1, v2, eps))
        return false;
    }
    return true;
  }
  
  public static boolean almostEquals(double[][] d1, double[][] d2, double eps) {
	    for (int i = 0; i < d1.length; i++) {
	    	for(int j=0;j<d1[0].length;j++){	    		
	    		double v1 = d1[i][j];
	    		double v2 = d2[i][j];
	    		if (!almostEquals(v1, v2, eps))
	    			return false;
	    	}
	    }
	    return true;
	  }
  public static double EPSILON = 0.000001;

  public static boolean almostEquals(double d1, double d2) {
    return almostEquals(d1, d2, EPSILON);
  }

  public static boolean almostEquals(double d1, double d2, double epsilon) {
    return Math.abs(d1 - d2) < epsilon;
  }

}