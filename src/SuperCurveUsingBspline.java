import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import javax.imageio.ImageIO;

public class SuperCurveUsingBspline {
	private static final String FILE_IN_PATH = "./data/data.txt";
	private static final String FILE_OUT_PATH = "./data/";

    BSplines bS;    
    
    public SuperCurveUsingBspline(){
        bS = new BSplines();
    }
    
    private double[][] mirrorW2d( double [][] s){

        double [][] s_mirror = new double[s.length+3][s[0].length+3];
        for(int i=0; i<s.length; i++){
            for(int j=0; j<s[0].length; j++){
                s_mirror[i+1][j+1] = s[i][j];
            }
        }
        
        /*########## mirror rows  1 and N-2 ################*/
        for(int j=0; j<s[0].length; j++){
            s_mirror[0][j+1] = s[1][j];
            s_mirror[s_mirror.length-2][j+1] = s[s.length-2][j];
        }
        
        /*########## mirror columns 1 and  N-2 ##############*/ 
        for(int i=0; i<s.length; i++){
            s_mirror[i+1][0] = s[i][1];
            s_mirror[i+1][s_mirror[0].length-2] = s[i][s[0].length-2];
        }   
        
        s_mirror[0][0] = s[1][1];
        s_mirror[0][s_mirror[0].length-2] = s[1][s[0].length-2];            
        s_mirror[s_mirror.length-2][0] = s[s.length-2][1];
        s_mirror[s_mirror.length-2 ][s_mirror[0].length-2] = 
                s[s.length-2][s[0].length-2];
                                        
        return s_mirror;
    }
    
    public double[][] cubicCoeff2d(double[][] s){
        DirectBsplFilter2d directFilter2d;
        directFilter2d = new DirectBsplFilter2d(s.length, s[0].length);
        double[][] coeffs = directFilter2d.filter(s);
        double[][] coeffs_mirror = mirrorW2d(coeffs);
        return coeffs_mirror;
    }
    public double MixingAndFusionOfTwoCurve(double[][] coeffs_mirror1,double[][] coeffs_mirror2, double row, double col){    
        int k = (int)Math.floor(row);   
        int l = (int)Math.floor(col);       
//        double interp_value =   
//        coeffs_mirror1[k+0][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+1)+ 
//        coeffs_mirror1[k+1][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+1)+
//        coeffs_mirror1[k+2][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+1)+
//        coeffs_mirror1[k+3][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+1)+
//                                                                    
//        coeffs_mirror1[k+0][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+0)+ 
//        coeffs_mirror1[k+1][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+0)+
//        coeffs_mirror1[k+2][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+0)+
//        coeffs_mirror1[k+3][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+0)+
//                                                                                                
//        coeffs_mirror1[k+0][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-1)+ 
//        coeffs_mirror1[k+1][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-1)+
//        coeffs_mirror1[k+2][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-1)+
//        coeffs_mirror1[k+3][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-1)+
//                                                
//        coeffs_mirror1[k+0][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-2)+ 
//        coeffs_mirror1[k+1][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-2)+
//        coeffs_mirror1[k+2][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-2)+
//        coeffs_mirror1[k+3][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-2);
        double interp_value =   
                coeffs_mirror1[k+0][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+1)+ 
                coeffs_mirror1[k+1][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+1)+
                coeffs_mirror1[k+2][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+1)+
                coeffs_mirror1[k+3][l+0]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+1)+
                                                                            
                coeffs_mirror1[k+0][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+0)+ 
                coeffs_mirror1[k+1][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+0)+
                coeffs_mirror1[k+2][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+0)+
                coeffs_mirror1[k+3][l+1]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+0)+
                                                                                                        
                coeffs_mirror1[k+0][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-1)+ 
                coeffs_mirror1[k+1][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-1)+
                coeffs_mirror1[k+2][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-1)+
                coeffs_mirror1[k+3][l+2]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-1)+
                                                        
                coeffs_mirror1[k+0][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-2)+ 
                coeffs_mirror1[k+1][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-2)+
                coeffs_mirror1[k+2][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-2)+
                coeffs_mirror1[k+3][l+3]*coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-2)+
                
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+1)+ 
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+1)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+1)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+1)+
                                                                            
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+0)+ 
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+0)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+0)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+0)+
                                                                                                        
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-1)+ 
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-1)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-1)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-1)+
                                                        
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-2)+ 
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-2)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-2)+
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-2);
        return interp_value;
    }
    
    public double[][] mixingOfCurveUsingBspline(double[][] s1,double[][] s2, int rate){
        double [][] coeffs_mirror1 = cubicCoeff2d(s1);
        double [][] coeffs_mirror2 = cubicCoeff2d(s2);
        int M1 = rate*s1.length - (rate-1);
        int N1 = rate*s1[0].length - (rate-1);
        double [][] s_interp1 = new double[M1][N1];
        for(int k=0; k<s_interp1.length; k++){
            for(int l=0; l<s_interp1[0].length; l++){
                s_interp1[k][l]=MixingAndFusionOfTwoCurve(coeffs_mirror1,coeffs_mirror2, k*(1.0/rate), l*(1.0/rate));
            }
        }                                                                                                       
        return s_interp1;                    
    }
    
	public static void main(String[] args){
    		//Write mixing txt file
            File out = new File(FILE_OUT_PATH +"out3.txt");
            double[][] pix_spec1=new double[3][72];
            double[][] pix_spec2=new double[3][72];
            double[][] MIXED_CURVE=new double[3][72];
    		try {
    			File file_1=new File("./data/avgLibrary.txt");
    			File file_2=new File("./data/avgLunarData.txt");
    		          Scanner sc1 = new Scanner(file_1);
    		          Scanner sc2 = new Scanner(file_2);
    		          sc1.useDelimiter(",");
    		          sc2.useDelimiter(",");
    		          for(int i=0;i<3;i++){
    		              for(int j=0;j<72;j++){
    		                 pix_spec1[i][j]=sc1.nextDouble();
    		                  //count++; 
    		                 pix_spec2[i][j]=sc2.nextDouble();
    		                 
    		                 System.out.print(pix_spec1[i][j]+"  ");
    		                 System.out.println(pix_spec2[i][j]);    		                 
    		              }
    		              System.out.println();
    		              sc1.nextLine();
    		              sc2.nextLine();
    		          }
    		          sc1.close();
    		          sc2.close();
    		          SuperCurveUsingBspline mixCurve = new SuperCurveUsingBspline();
    		          MIXED_CURVE=mixCurve.mixingOfCurveUsingBspline(pix_spec1,pix_spec2,3);
    		          for(int i=0;i<3;i++){
    		              for(int j=0;j<72;j++){    		                 
    		                 System.out.print(MIXED_CURVE[i][j]+"  ");    		                 
    		              }
    		              System.out.println();
    		          }
//    			PrintWriter pw = new PrintWriter(new FileWriter(out));
//    			pw.close();
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
        }
    
}