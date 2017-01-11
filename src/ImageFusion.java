import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

public class ImageFusion {
    BSplines bS;    
    
    public ImageFusion(){
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
	public static BufferedImage[] FuseTwoImages(BufferedImage[] hysi,BufferedImage tmc){
		BufferedImage[] finalImage=new BufferedImage[65];
		for(int k=1;k<=64;k++){
			int height=hysi[k].getHeight();
			int width=hysi[k].getWidth();
			double red1[][]=new double[height][width];
			double green1[][]=new double[height][width];
			double blue1[][]=new double[height][width];
			double red2[][]=new double[height][width];
			double green2[][]=new double[height][width];
			double blue2[][]=new double[height][width];
			for(int i=0;i<height;i++){
					for(int j=0;j<width;j++){   
						Color c=new Color(hysi[k].getRGB(j, i));
						red1[i][j]=c.getRed();
						blue1[i][j]=c.getBlue();
						green1[i][j]=c.getGreen();
					}
				}
			 for(int i=0;i<hysi[1].getHeight();i++){
					for(int j=0;j<width;j++){   
						Color c=new Color(tmc.getRGB(j, i));
						red2[i][j]=c.getRed();
						blue2[i][j]=c.getBlue();
						green2[i][j]=c.getGreen();
					}
				}
				ImageFusion mixCurve = new ImageFusion();
				ImageFusion mixCurve1 = new ImageFusion();
				ImageFusion mixCurve2 = new ImageFusion();
				double[][] MIXED_CURVE1=mixCurve.mixingOfCurveUsingBspline(red1,red2,1);
				double[][] MIXED_CURVE2=mixCurve1.mixingOfCurveUsingBspline(green1,green2,1);
				double[][] MIXED_CURVE3=mixCurve2.mixingOfCurveUsingBspline(blue1,blue2,1);
				finalImage[k]=ColorToImages(MIXED_CURVE1,MIXED_CURVE2,MIXED_CURVE3);
				System.out.println("complete finish"+k);
		}
		return finalImage;
	}
	public static BufferedImage ColorToImages(double[][] red,double[][] green,double[][] blue){
		BufferedImage check=new BufferedImage(red[0].length,red.length,BufferedImage.TYPE_INT_RGB);
		 for (int y = 0; y < red.length; y++) {
		     for (int x = 0; x < red[0].length; x++) {
		        int rgb = (int)red[y][x];
		        rgb = (rgb << 8) + (int)green[y][x]; 
		        rgb = (rgb << 8) + (int)blue[y][x];
		        check.setRGB(x, y, rgb);
		     }
		  }
		return check;
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
        double interp_value =   
                (coeffs_mirror1[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+1)+ 
                coeffs_mirror1[k+1][l+0]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+1)+
                coeffs_mirror1[k+2][l+0]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+1)+
                coeffs_mirror1[k+3][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+1)+
                                                                            
                coeffs_mirror1[k+0][l+1]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+0)+ 
                coeffs_mirror1[k+1][l+1]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l+0)+
                coeffs_mirror1[k+2][l+1]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l+0)+
                coeffs_mirror1[k+3][l+1]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l+0)+
                                                                                                        
                coeffs_mirror1[k+0][l+2]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-1)+ 
                coeffs_mirror1[k+1][l+2]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-1)+
                coeffs_mirror1[k+2][l+2]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-1)+
                coeffs_mirror1[k+3][l+2]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-1)+
                                                        
                coeffs_mirror1[k+0][l+3]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l-2)+ 
                coeffs_mirror1[k+1][l+3]*bS.bspline(3,row-k+0)*bS.bspline(3,col-l-2)+
                coeffs_mirror1[k+2][l+3]*bS.bspline(3,row-k-1)*bS.bspline(3,col-l-2)+
                coeffs_mirror1[k+3][l+3]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-2))+
                
                (coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k+1)*bS.bspline(3,col-l+1)+ 
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
                coeffs_mirror2[k+0][l+0]*bS.bspline(3,row-k-2)*bS.bspline(3,col-l-2))/2;
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
    public static BufferedImage blendHysi(BufferedImage[] bi, double weight) {
//		if (bi == null) throw new NullPointerException("bi1 is null");

//		if (bi2 == null) throw new NullPointerException("bi2 is null");

		int width = bi[1].getWidth();
//		if (width != bi2.getWidth()) throw new IllegalArgumentException("widths not equal");

		int height = bi[1].getHeight();
//		if (height != bi2.getHeight()) throw new IllegalArgumentException("heights not equal");

		BufferedImage bi3 = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		int rgbim[][] = new int[65][width];
		int rgbim65[] = new int[width];

		for (int row = 0; row < height; row++) {
			for(int l=1;l<=64;l++){
				bi[l].getRGB(0, row, width, 1, rgbim[l], 0, width);				
			}
			int[] r=new int[65];
			int[] g=new int[65];
			int[] b=new int[65];
			for (int col = 0; col < width; col++) {
				int rgb[]=new int[65];
				for(int k=1;k<=64;k++){					
					rgb[k] =rgbim[k][col];
					r[k] = (rgb[k] >> 16) & 255;
					g[k] = (rgb[k] >> 8) & 255;
					b[k] = rgb[k] & 255;
				}
				int r65=0,g65=0,b65=0;
				for(int m=1;m<=64;m++){					
					r65 += (int)(r[m] * weight);
					g65 += (int)(g[m] * weight);
					b65 += (int)(b[m] * weight);
				}
				rgbim65[col] = (r65 << 16) | (g65 << 8) | b65;
			}

			bi3.setRGB(0, row, width, 1, rgbim65, 0, width);
		}

		return bi3;
	}

    public BufferedImage doubleArrayToImage(double[][] array){     
        BufferedImage image = new BufferedImage(array[0].length,array.length,BufferedImage.TYPE_INT_RGB);  
        for( int y = 0; y < array.length; y++ ){
            for( int x = 0; x < array[0].length; x++ ){ 
                 int value = (int)array[y][x] << 16|(int)array[y][x] << 8|(int)array[y][x];
                 image.setRGB(x, y, value);                                                   
            }
        }       
        return image;
    }
    
	public static void main(String[] args){
		BufferedImage image=null;
        BufferedImage[] fuse=new BufferedImage[65];
        BufferedImage[] Hysi_Img=new BufferedImage[65];
        BufferedImage TMC_IMG=null;
        try {
        	TMC_IMG=ImageIO.read(new File("./data/tmc2561.png"));
        	for(int j=1;j<=64;j++){
        		Hysi_Img[j]=ImageIO.read(new File("./data/himg/h"+j+".png"));
        		System.out.println("take"+j);
	        }
            fuse=ImageFusion.FuseTwoImages(Hysi_Img, TMC_IMG);
            ImageView imageView4 = new ImageView();
    		ImageView imageView3 = new ImageView();
    		ImageView imageView1 = new ImageView();
    		ImageView imageView5 = new ImageView();
    		imageView4.drawImage(fuse[2]);
    		imageView3.drawImage(fuse[3]);
    		imageView1.drawImage(fuse[6]);
    		image=ImageFusion.blendHysi(fuse,0.015);
    		imageView5.drawImage(image);
    		ImageIO.write(image, "PNG", new File("./data/Fused.png"));
		} catch (IOException e) {
				e.printStackTrace();
		}
    }
    public class BSplines {

        public double bspline(int degree, double x){
            double betta;
            double t;
            betta = 0;
            if(degree == 0){                
                if ((x > -0.5) && (x < 0.5)){
                    betta = 1.0;
                }
                else if( Math.abs(x) == 0.5){
                    betta = 0.5;
                }
                else if( Math.abs(x) > 0.5){
                    betta = 0.0;
                }           
            }
            else if( degree == 1){
                if ((x<=-1) || (x>=1)){ 
                    betta = 0.0;                        
                }
                else if ((x>-1) && (x<0)){
                    betta = x+1;
                }
                else if ((x>0) && (x<1)){
                    betta = -x+1;
                }
                else if( x==0){
                    betta = 1.0;
                }                                   
            }       
            else if (degree == 2 ){     
                t = 1.5;
                if ((x <= (0-t)) || (x >= (3-t))){
                    betta = 0.0;
                }
                else if ((x >= (0-t)) && (x< (1-t))) {
                    betta = ((x+t)*(x+t))/2.0;
                }
                else if ((x >= (1-t)) && (x< (2-t))) {
                    betta = ((x+t)*(x+t)-3.0*(x-1+t)*(x-1+t))/2.0;
                }
                else if ((x >= (2-t)) && (x< (3-t))) {
                    betta = ((x+t)*(x+t) - 3.0*(x-1+t)*(x-1+t) + 
                            3.0*(x-2+t)*(x-2+t))/2.0;
                }
            }
            else if (degree == 3 ){ 
                if ((Math.abs(x)>=0) && (Math.abs(x)< 1)) {
                    betta = 2.0/3.0 - Math.abs(x)*Math.abs(x) + 
                    (Math.abs(x)*Math.abs(x)*Math.abs(x))/2.0;
                }
                else if ((Math.abs(x)>=1) && (Math.abs(x)< 2)) {
                    betta = ((2-Math.abs(x))*(2-Math.abs(x))*
                            (2-Math.abs(x)))/6.0;
                }
                else if (Math.abs(x) >=2) {
                    betta = 0.0;
                }
            }
            return betta;
        }
    }
    public class DirectBsplFilter1d {

        private int N;
        private double z1;
        private double cplus [];
        private double cminus [];
        private int k;
        private double sum0;
        
        public DirectBsplFilter1d(int N){
            this.N = N;
            z1 = -2.0 +  Math.sqrt(3.0);
            cplus = new double[N];
            cminus = new double[N];
        }
        
        public void reset(){
            for(k=0; k<N; k++){
                cplus[k] = 0.0;
                cminus[k] = 0.0;
            }
        }
        public double [] filter(double s[]){
            sum0 = 0.0;
            for(k = 0; k < N; k++){
                sum0 = sum0 + 6.0*s[k]*Math.pow(z1, k);
            }
            cplus[0] = sum0;
            for(k = 1; k < N; k++){
                cplus[k] = 6.0*s[k] + z1*cplus[k-1];
            }
            cminus[N-1] = (z1/(z1*z1-1.0))*
                    (cplus[N-1] + z1*cplus[N-2]);
            for(k = N-2; k >= 0; k--){
                cminus[k] = z1*(cminus[k+1]-cplus[k]);
            }
            return cminus;
        }
    }
    public class DirectBsplFilter2d {

        private DirectBsplFilter1d directFilter1dX; 
        private DirectBsplFilter1d directFilter1dY;         
        
        public DirectBsplFilter2d(int M, int N){
            directFilter1dX = new DirectBsplFilter1d(M);    
            directFilter1dY = new DirectBsplFilter1d(N);    
        }
        
        public double[][] filter(double [][] img){  
            double [] row = new double[img[0].length];
            double [] col = new double[img.length];
            double [] filt_row = new double[img[0].length];
            double [] filt_col = new double[img.length];
            double [][] coeffs = new double[img.length][img[0].length];
            
            /*#################### filtrations along y ##################*/
            for(int i=0; i<img.length; i++){
                for(int j=0; j<img[0].length; j++){
                    row[j] = img[i][j];
                }
                filt_row = directFilter1dY.filter(row);
                for(int j=0; j<img[0].length; j++){
                    coeffs[i][j] = filt_row[j];
                }
                directFilter1dY.reset();
            }
            /*#################### filtrations along y ##################*/
                        
            /*#################### filtrations along x ##################*/
            for(int j=0; j<img[0].length; j++){
                for(int i=0; i<img.length; i++){
                    col[i] = coeffs[i][j];
                }
                filt_col = directFilter1dX.filter(col);
                for(int i=0; i<img.length; i++){
                    coeffs[i][j] = filt_col[i];
                }
                directFilter1dX.reset();
            }
            /*#################### filtrations along x ##################*/             
            return coeffs;                                                  
        }

    }    
}