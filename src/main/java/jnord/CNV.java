package jnord;

import java.util.Arrays;

public class CNV extends RODImpl {
    String cnvType;
    Gene gene;
    public CNV(String chr, int start, int end, String cnvType, Gene g){
        super(chr, start, end);
        this.cnvType = cnvType;
        this.gene = g;
    }
    public String getType(){
        return this.cnvType;
    }
    public void setType(String t){
        this.cnvType = t;
    }
    public Gene getGene(){
        return this.gene;
    }
    public double[] getRatioNormalized(String sample, int cnvStart, int cnvEnd){
        double[] buf = new double[cnvEnd - cnvStart + 1];
        int len = 0;
        for(Region region: this.gene.getRegions()){
            // System.out.println("about region from " + region.getStart() + " to " + region.getEnd());
            double[] normalized = region.getRatioNormalized(sample);
            for(int pos = region.getStart(); pos < region.getEnd(); pos++){
                if(pos >= cnvStart && pos < cnvEnd){
                    // System.out.println("cnv len[" + buf.length+ "], region len[" + normalized.length + "], pos["+ pos + "], index[" + (pos - region.getStart()) + "] " + sample + " " + gene.getName());
                    buf[len] = normalized[pos - region.getStart()];
                    len++;
                }
            }
        }
        double[] result = Arrays.copyOfRange(buf, 0, len);
        return result;
    }
    public double[] getMedian(int cnvStart, int cnvEnd){
        double[] buf = new double[cnvEnd - cnvStart + 1];
        int len = 0;
        for(Region region: this.gene.getRegions()){
            // System.out.println("about region from " + region.getStart() + " to " + region.getEnd());
            double[] median = region.getMedian();
            for(int pos = region.getStart(); pos < region.getEnd(); pos++){
                if(pos >= cnvStart && pos < cnvEnd){
                    buf[len] = median[pos - region.getStart()];
                    len++;
                }
            }
        }
        double[] result = Arrays.copyOfRange(buf, 0, len);
        return result;

    }
}
