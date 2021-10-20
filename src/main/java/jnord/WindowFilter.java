package jnord;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

public class WindowFilter {
    double passThreshold = 180;
    double lengthThreshold = 200;
    public WindowFilter(double pass, double length){
        this.passThreshold = pass;
        this.lengthThreshold = length;
    }
    public WindowFilter(){
    }
    public boolean accept(double[] buf, double[] median, Detector detector){
        int passCount = 0;
        for(int i = 0; i<buf.length; i++){
            if(detector.check(buf[i], median[i])>1){
                passCount++;
            }
        }
        if(buf.length > lengthThreshold && passCount > passThreshold){
            return true;
        }
        return false;
    }
}
