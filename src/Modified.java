import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

public class Modified {

	public static void main(String[] args) {
		try {
			ImageView imageView = new ImageView();
			ImageView imageView2 = new ImageView();
			ImageView imageView3 = new ImageView();
			BufferedImage bi1 = ImageIO.read(new File("./data/take1.png"));
			BufferedImage bi2 = ImageIO.read(new File("./data/take2.png"));
		 	int height=bi1.getHeight();
			  int width=bi2.getWidth();	
			  double red1[][]=new double[height][width];
			 double green1[][]=new double[height][width];
			 double blue1[][]=new double[height][width];
			 for(int i=0;i<bi1.getHeight();i++)
				{
					for(int j=0;j<bi1.getWidth();j++)
					{   
						Color c=new Color(bi1.getRGB(j, i));
						red1[i][j]=c.getRed();
						blue1[i][j]=c.getBlue();
						green1[i][j]=c.getGreen();
					}
				}
			 double red2[][]=new double[height][width];
			 double green2[][]=new double[height][width];
			 double blue2[][]=new double[height][width];
			 for(int i=0;i<bi1.getHeight();i++)
				{
					for(int j=0;j<bi1.getWidth();j++)
					{   
						Color c=new Color(bi2.getRGB(j, i));
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
			 BufferedImage check=new BufferedImage(red1[0].length,red1.length,BufferedImage.TYPE_INT_RGB);
			 for (int y = 0; y < MIXED_CURVE1.length; y++) {
			     for (int x = 0; x < MIXED_CURVE1[0].length; x++) {
			        int rgb = (int)MIXED_CURVE1[y][x];
			        rgb = (rgb << 8) + (int)MIXED_CURVE2[y][x]; 
			        rgb = (rgb << 8) + (int)MIXED_CURVE3[y][x];
			        check.setRGB(x, y, rgb);
			     }
			  }
			imageView.drawImage(bi1);
			imageView2.drawImage(bi2);
			imageView3.drawImage(check);
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

}
