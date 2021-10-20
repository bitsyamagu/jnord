package jnord;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ConcurrentHashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;

public class TableWriter {
    ArrayList<String> samples;
    ArrayList<Gene> genes;
    List<String> geneNames;
    String[] header = {
        "#chr",
        "pos",
        "gene",
        "median",
        "sd",
        "sn",
        "normalized_sn"
    };
    public TableWriter(List<String> names, ArrayList<Gene> genes, ArrayList<String> samples){
        this.genes = genes;
        this.samples = samples;
        this.geneNames = names;
    }
    public static String join(String[] buf){
        StringBuilder sb = new StringBuilder();
        sb.append(buf[0]);
        for(int i = 1; i<buf.length; i++){
            sb.append("\t");
            sb.append(buf[i]);
        }
        return sb.toString();
    }
    public void write() throws IOException{
        HashMap<String, Gene> map = new HashMap<>();
        for(Gene gene: this.genes){
            map.put(gene.getName(), gene);
        }
        for(String geneName: this.geneNames){
            Gene gene = map.get(geneName);
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(gene.getName() + ".bed")));
            out.print(join(header));
            for(String sample: this.samples){
                out.print("\t");
                out.print(sample + ".depth");
                out.print("\t");
                out.print(sample + ".normzlied");
                out.print("\t");
                out.print(sample + ".normalized_ratio");
            }
            out.println("");
            for(Region region: gene.getRegions()){
                int start = region.getStart();
                int end = region.getEnd();

                double[] median = region.median;
                double[] sd = region.sd;
                double[] sn = region.sn;
                double[] normalizedSN = region.normalizedSN;

                ConcurrentHashMap<String, double[]> depth = region.depthMap;
                ConcurrentHashMap<String, double[]> normalized = region.normalized;
                ConcurrentHashMap<String, double[]> ratioNormalized = region.ratioNormalized;
                for(int pos = start; pos < end; pos++){
                    out.print(region.getChromosome());
                    out.print("\t");
                    out.print(pos);
                    out.print("\t");
                    out.print(gene.getName());
                    out.print("\t");
                    out.print(median[pos - start]);
                    out.print("\t");
                    out.print(sd[pos - start]);
                    out.print("\t");
                    out.print(sn[pos - start]);
                    out.print("\t");
                    out.print(normalizedSN[pos - start]);
                    for(String sample: this.samples){
                        out.print("\t");
                        out.print(depth.get(sample)[pos - start]);
                        out.print("\t");
                        out.print(normalized.get(sample)[pos - start]);
                        out.print("\t");
                        out.print(ratioNormalized.get(sample)[pos - start]);
                    }
                    out.println("");
                }
            }
            out.close();
        }
    }

}
