/**
 * see 
 *  https://rdrr.io/bioc/affy/man/normalize.invariantset.html
 *  https://rdrr.io/bioc/affy/src/R/normalize.invariantset.R
 */
package jnord;

import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

public class InvariantSetNormalization {
    PolynomialSplineFunction spline;
    PolynomialSplineFunction extrapolatorSpline;
    PolynomialFunction firstSpline = null;
    PolynomialFunction lastSpline = null;
    // PolynomialFunction[] curve = null;
    double[] range = null;
    boolean debug = false;

    private void hist(double[] buf){
        int[] hist = new int[11];
        for(int i = 0; i<buf.length; i++){
            int index = (int)(buf[i]/0.1);
            hist[index]++;
        }
   //   for(int i = 0; i<hist.length; i++){
   //       System.out.println(i + " " + hist[i]);
   //   }
    }
    public void setDebug(boolean b){
        this.debug = b;
    }
    double prd_td_min = 0.003;
    double prd_td_max = 0.007;
    public void buildSpline(double[] data, double[] median){
        NaturalRanking ranking = new NaturalRanking();
        int np = data.length;
        double[] rank_ref = ranking.rank(median);
        double[] rank_array = ranking.rank(data);

        // Init --
        // adjusted threshold things
        double prd_td_min_adj = prd_td_min * 10; 
        double prd_td_max_adj = prd_td_max * 10;

        // index all the PM probes as being in the invariant set
        boolean[] iset = new boolean[np];
        for(int i = 0; i<np; i++){
            iset[i] = true;
        }
        // number of probes in the invariant set
        int ns = np;
        // number of probes previously in the invariant set
        int nsOld = ns + 50 + 1;
        // End of Init --

        // iterate while the number of genes in the invariant set (ns) still varies...
        while((nsOld - ns) > 50){
     //     System.err.println("(nsOld - ns) > 50 : (" + nsOld + " - " + ns + ") > 50");
            // average intensity rank for the probe intensities
            double[] air = divide( add(filt(rank_ref, iset), filt(rank_array, iset)), 2*ns);
            double[] prd = divide(
                abs(
                    sub( filt(rank_ref, iset), filt(rank_array, iset))
                    ), (double)ns);
            hist(prd);
            double[] threshold = add(multi(air, (prd_td_max_adj - prd_td_min)), prd_td_min_adj);
            iset  = lessThanAnd(prd, iset, threshold);

            nsOld = ns;
            ns = sum(iset);

            // update adjusted threshold parameters
            if(prd_td_min_adj > prd_td_min){ 
                prd_td_min_adj *= 0.9;
                prd_td_max_adj *= 0.9;
            }
     //     System.err.println("(nsOld - ns) > 50 : (" + nsOld + " - " + ns + ") > 50");
        }
        /*
          ncurve = smooth.spline(filt(median, iset), filt(data, iset));
        */
        LoessInterpolator intp = new LoessInterpolator();
        // SplineInterpolator intp = new SplineInterpolator();
        // double[] repr_median = filt(median, iset);
        // double[] repr_data = filt(data, iset);

        double[] hoge1 = filt(data, iset);
        double[] hoge2 = filt(median, iset);
        double[][] plot = sort(hoge1, hoge2);
  //    for(int i = 0; i<plot[0].length; i++){
  //          System.out.println(plot[0][i] + " " + plot[1][i]);
  //    }

        this.spline = intp.interpolate(plot[1], plot[0]);
        this.range = spline.getKnots();
        System.err.println("-- Spline interpolation finished --");
        System.err.println("From: " + this.range[0]);
        System.err.println("To  : " + this.range[this.range.length-1]);
        System.err.println("--");

        PolynomialFunction[] curve = null;
        curve = this.spline.getPolynomials();
        // PolynomialFunction:
        this.firstSpline = curve[2];
        this.lastSpline = curve[curve.length-4];
        // 
    //  PolynomialFunction f1 = curve[curve.length-3];
    //  System.out.println("degree: " + f1.degree() + " " + java.util.Arrays.toString(f1.getCoefficients()));
    //  PolynomialFunction f2 = curve[curve.length-2];
    //  System.out.println("degree: " + f2.degree() + " " + java.util.Arrays.toString(f2.getCoefficients()));
    //  PolynomialFunction f3 = curve[curve.length-1];
    //  System.out.println("degree: " + f3.degree() + " " + java.util.Arrays.toString(f3.getCoefficients()));
    //  // for global
    //  double[] global_1 = new double[3];
    //  double[] global_2 = new double[3];
    //  global_1[0] = plot[0][0];
    //  global_1[1] = plot[0][plot[0].length/2];
    //  global_1[2] = plot[0][plot[0].length-1];
    //  global_2[0] = plot[1][0];
    //  global_2[1] = plot[1][plot[1].length/2];
    //  global_2[2] = plot[1][plot[1].length-1];
    //  SplineInterpolator intp2 = new SplineInterpolator();
    //  this.extrapolatorSpline = intp2.interpolate(global_1, global_2);
    //  PolynomialFunction[] funcs = extrapolatorSpline.getPolynomials();
    //  this.firstSpline = funcs[0];
    //  this.lastSpline = funcs[funcs.length-1];
    //  this.range = this.extrapolatorSpline.getKnots();
        // this.lastSpline = funcs[funcs.length-1];
    }
    public PolynomialSplineFunction getCurve(){
        // return this.extrapolatorSpline;
        return this.spline;
    }
    class Point implements Comparable<Point> {
        double x;
        double y;
        Point(double x_, double y_){
            x = x_;
            y = y_;
        }
        @Override
        public int compareTo(final Point another){
            int c = Double.compare(y, another.y);
            if(c==0){
                return Double.compare(x, another.x);
            }
            return c;
        }
    }
    // sort and collapse points which has same X
    public double[][] sort(double[] buf1, double[] buf2){
        ArrayList<Point> list = new ArrayList<>();
        for(int i = 0; i<buf1.length; i++){
            list.add(new Point(buf1[i], buf2[i]));
        }
        Collections.sort(list);

        ArrayList<Point> list2 = new ArrayList<>();
        Point last = list.get(0);
        for(int i = 0; i<buf1.length; i++){
            Point p = list.get(i);
            if(last.y == p.y){
                continue;
            }else {
                list2.add(p);
            }
            last = p;
        }
        double[][] ret = new double[2][];
        ret[0] = new double[list2.size()];
        ret[1] = new double[list2.size()];
        for(int i = 0; i<list2.size(); i++){
            ret[0][i] = list2.get(i).x;
            ret[1][i] = list2.get(i).y;
        }
        return ret;
    }
    public static boolean[] lessThanAnd(double[] buf, boolean[] flags, double[] threshold){
        boolean[] ret = Arrays.copyOf(flags, flags.length);
        int index = 0;
        for(int i = 0; i<flags.length; i++){
            boolean flag = flags[i];
            if(flag){
                if(buf[index] < threshold[index]){
                    ret[i] = true;
                }else {
                    ret[i] = false;
                }
                index++;
            }
        }
        return ret;
    }
    public static int sum(boolean[] flags){
        int c = 0;
        for(boolean b: flags){
            if(b){ c++; }
        }
        return c;
    }
    public static double[] abs(double[] buf){
        double[] ret = new double[buf.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = Math.abs(buf[i]);
        }
        return ret;
    }
    public static double[] multi(double[] buf, double c){
        double[] ret = new double[buf.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = buf[i] * c;
        }
        return ret;
    }
    public static double[] add(double[] buf, double c){
        double[] ret = new double[buf.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = buf[i] + c;
        }
        return ret;
    }
    public static double[] divide(double[] buf, double denom){
        double[] ret = new double[buf.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = buf[i] / denom;
        }
        return ret;
    }
    public static double[] add(double[] buf1, double[] buf2){
        double[] ret = new double[buf1.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = buf1[i] + buf2[i];
        }
        return ret;
    }
    public static double[] sub(double[] buf1, double[] buf2){
        double[] ret = new double[buf1.length];
        for(int i = 0; i<ret.length; i++){
            ret[i] = buf1[i] - buf2[i];
        }
        return ret;
    }
    public synchronized static double[] filt(double[] buf, boolean[] flags){
        int count = 0;
        for(boolean b: flags){
            if(b){
                count++;
            }
        }
        double[] ret = new double[count];
        count = 0;
        for(int i = 0; i<buf.length; i++){
            if(flags[i]){
                ret[count] = buf[i];
                count++;
            }
        }
        return ret;
    }
    public void init(double[] rawIntensities, double[] rawReference){
        double[] logRaw = new double[rawIntensities.length];
        double[] logRef = new double[rawReference.length];
        // convert and fix -Infinity
        for(int i = 0; i<logRaw.length; i++){
            logRaw[i] = Math.log(rawIntensities[i]);
            if(Double.isInfinite(logRaw[i])){ logRaw[i] = 0; };
        }
        for(int i = 0; i<logRef.length; i++){
            logRef[i] = Math.log(rawReference[i]);
            if(Double.isInfinite(logRef[i])){ logRef[i] = 0; };
        }
        buildSpline(logRaw, logRef);
    }
    double predict(double x, double[] coefficients){
        return coefficients[3]*x*x*x + coefficients[2]*x*x + coefficients[1]*x + coefficients[0];
    }
    public double normalize(double raw, double ref){
        double logRaw = Math.log(raw);
        double logRef = Math.log(ref);
        
        /// if(logRef <= 0 || logRef == Double.NEGATIVE_INFINITY){
        if(logRef == Double.NEGATIVE_INFINITY){
            // logRef = 0;
            return 1.0;
        }
        if(logRaw == Double.NEGATIVE_INFINITY){
            // logRaw = logRef; // treat it as median(control)
            // logRaw = 0;
            return 1.0;
            // return Math.exp(firstSpline.value(logRef));
        }
        if(logRef < this.range[0]){
            // System.err.println("too small intensity: " + logRef + " " + firstSpline.value(logRef));
            // System.out.println("logRaw:logRef:firstSpline.value(logRef) "+ logRaw + " " + logRef + " "+ firstSpline.value(logRef));
            // return Math.exp(logRaw + logRef - firstSpline.value(logRef));
            return Math.exp(logRaw);
        }else if(logRef > this.range[this.range.length-1]){
            // System.out.println("logRaw:logRef:lastSpline.value(logRef) "+ logRaw + " " + logRef + " "+ lastSpline.value(logRef));
            return Math.exp(logRaw + logRef - lastSpline.value(logRef));
        }
        try {
            if(debug){
                System.out.println("type 3: " + spline.value(logRef));
            }
            // spline.value(logRef);
        }catch(Exception e){
            e.printStackTrace();
        }
        return Math.exp(logRaw + logRef - spline.value(logRef));
    }
}
