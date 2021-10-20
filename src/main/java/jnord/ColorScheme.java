package jnord;

import java.awt.Color;

public class ColorScheme {
    static int saturation = 255;
    public static void setSaturation(int s){
        saturation = s;
    }
    static double lastH = 0.0;
    public static Color randomColor(){
        double golden_ratio_conjugate = 0.618033988749895;
        double h = (Math.random() + golden_ratio_conjugate) % 1 *360;
        while(Math.abs(h - lastH) < 80.0){
            h = (Math.random() + golden_ratio_conjugate) % 1 *360;
        }
        // System.err.println("color diff: " + h + " "+ (h - lastH));
        lastH = h;
        return hsvToRgb(h, 50, 95);
        // return "rgb("+rgb[0]+","+rgb[1]+","+rgb[2]+")";
    }

    /**
    * Converts an HSV color value to RGB. Conversion formula
    * adapted from http://en.wikipedia.org/wiki/HSL_and_HSV.
    * Assumes h is contained in the set [0, 360] and
    * s and l are contained in the set [0, 100] and
    * returns r, g, and b in the set [0, 255].
    *
    * param   Number  h       The hue
    * param   Number  s       The saturation
    * param   Number  v       The value
    * return  Array           The RGB representation
    */
    private static Color hsvToRgb(double h, double s, double v){
        double chroma = s * v / 10000;
        double min = v / 100 - chroma;
        double hdash = h / 60;
        double x = chroma * (1 - Math.abs(hdash % 2 - 1));
        double r = 0, g = 0, b = 0;

        
        if(hdash < 1){
          r = chroma;
          g = x;
        }else if (hdash < 2){
          r = x;
          g = chroma;
        }else if (hdash < 3){
          g = chroma;
          b = x;
        }else if (hdash < 4){
          g = x;
          b = chroma;
        }else if (hdash < 5){
          r = x;
          b = chroma;
        }else if (hdash <= 6){
          r = chroma;
          b = x;
        }
        
        r += min;
        g += min;
        b += min;
        
        /*
        System.err.println("" + r + " " + g + " " +  b);
        System.err.println("" + 
            (int)Math.round(r * 255) + " " +
            (int)Math.round(g * 255) + " " +
            (int)Math.round(b * 255));
        */
        return new Color(
            (int)Math.round(r * saturation),
            (int)Math.round(g * saturation),
            (int)Math.round(b * saturation)
        );
    }
    public static void main(String[] argv){
        System.setProperty("java.awt.headless", "true");
        for(int i = 0; i<10; i++){
            System.out.println(randomColor().toString());
        }
    }
}
