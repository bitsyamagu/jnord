package jnord;

import java.util.concurrent.ConcurrentHashMap;
import java.util.Collection;

public class Region extends RODImpl {
    public String gene;
    public String nm;
    public double[] median;
    public double[] gcContent;
    public double[] sd;
    public double[] sn;
    public double[] normalizedSN;
    ConcurrentHashMap<String, double[]> depthMap = new ConcurrentHashMap<>();
    ConcurrentHashMap<String, double[]> normalized = new ConcurrentHashMap<>();
    ConcurrentHashMap<String, double[]> ratioNormalized = new ConcurrentHashMap<>();
    public Region(String c, int s, int e, String g, String n){
        super(c, s, e);
        gene = g;
        nm = n;
        median = new double[e-s + 1];
    }
    public void addDepth(String sampleName, double[] buf){
        // synchronized(depthMap){
            depthMap.put(sampleName, buf);
        // }
    }
    public double[] getDepth(String sampleName){
        return this.depthMap.get(sampleName);
    }
    public Collection<double[]> getDepths(){
        return depthMap.values();
    }
    public void addNormalized(String sampleName, double[] buf){
        normalized.put(sampleName, buf);
    }
    public void addRatioNormalized(String sampleName, double[] buf){
        ratioNormalized.put(sampleName, buf);
    }
    public double[] getNormalizedBySample(String sample){
        return this.normalized.get(sample);
    }
    public Collection<double[]> getNormalized(){
        return normalized.values();
    }
    public double[] getRatioNormalized(String sample){
        return this.ratioNormalized.get(sample);
    }
    public void setMedian(double[] buf){
        this.median = buf;
        // dump();
    }
    public double[] getMedian(){
        return this.median;
    }
    public void setSD(double[] buf){
        this.sd = buf;
    }
    public double[] getSD(){
        return this.sd;
    }
    public void setSN(double[] buf){
        this.sn = buf;
    }
    public double[] getSN(){
        return this.sn;
    }
    public void setNormalizedSN(double[] buf){
        this.normalizedSN = buf;
    }
    public double[] getNormalizedSN(){
        return this.normalizedSN;
    }
    public void setGCContent(double[] buf){
        this.gcContent = buf;
    }
    public double[] getGCContent(){
        return this.gcContent;
    }
    public void dump(){
        for(int i = 0; i<this.median.length; i++){
            System.err.print(this.gene + " "+ getChromosome() + " " + (getStart() + i) + " " + median[i]);
            for(double[] buf: depthMap.values()){
                System.err.print("\t" + buf[i]);
            }
            System.err.println("");
        }
    }
}
