package jnord;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

public class MedianFilter {
    double lossThreshold = 0.7;
    double gainThreshold = 1.3;
    public MedianFilter(double loss, double gain){
        this.lossThreshold = loss;
        this.gainThreshold = gain;
    }
    public MedianFilter(){
    }
    public boolean accept(double[] buf){
        Percentile perc = new Percentile();
        double median = perc.evaluate(buf, 50.0f);
        if(median <= lossThreshold || median >= gainThreshold){
            return true;
        }
        return false;
    }
}
