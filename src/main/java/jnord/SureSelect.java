package jnord;

import java.util.List;
import java.util.HashSet;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Collections;

/**
 * This object extracts regions for target genes
 */

public class SureSelect {
    HashSet<String> targetGenes = new HashSet<String>();
    List<Region> regions = new ArrayList<Region>();
    boolean trimChromosomePrefix = true;
    /*
    public class Region extends RODImpl {
        public String gene;
        public String nm;
        public Region(String c, int s, int e, String g, String n){
            super(c, s, e);
            gene = g;
            nm = n;
        }
    };*/
    public SureSelect(){
    }
    public SureSelect(File file, List<String> targets, boolean trimChr) throws FileNotFoundException, IOException{
        this.trimChromosomePrefix = trimChr;
        this.targetGenes.addAll(targets);
        this.regions = parse(file);
    }
    public List<Region> getRegions(){
        return this.regions;
    }
    public List<Region> getRegions(String geneName){
        ArrayList<Region> list = new ArrayList<Region>();
        for(Region r: this.regions){
            // System.out.println("in getRegion(): "+ list.size() );
            if(r.gene.equals(geneName)){
                list.add(r);
                // System.out.println("//" + geneName + " " + r.getStart() + " " + r.getEnd());
            }
        }
        return list;
    }
    public List<Region> parse(File file) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(file));
        String raw = null;
        Pattern refPattern = Pattern.compile("^ref\\|(\\S+)");
        Pattern nmPattern = Pattern.compile("^ref\\|(NM_[0-9]+)");
        while(null != (raw = br.readLine())){
            if(!raw.startsWith("chr")){
                continue;
            }
            String[] line = raw.split("\t", 4);
            String[] info = line[3].split(",");
            Matcher nmmat = nmPattern.matcher(line[3]);
            String nm = "";
            if(nmmat.matches()){
                nm = nmmat.group(1);
            }
            for(String inf: info){
                Matcher m = refPattern.matcher(inf);
                if(m.matches()){
                    String gene = m.group(1); // 0 is whole of hit
                    if(this.targetGenes.contains(gene)){
                        // System.out.println(raw);
                        if(trimChromosomePrefix){
                            if(line[0].startsWith("chr") || line[0].startsWith("Chr")){
                                line[0] = line[0].substring(3);
                            }
                        }
                        regions.add(
                            new Region(line[0], Integer.parseInt(line[1]), Integer.parseInt(line[2]), gene, nm)
                        );
                    }
                }
            }
        }
        br.close();
        merge();
    //  for(Region r: regions){
    //          System.out.println("SureSelect: " + r.gene + "\t" + r.chr + "\t" + r.start + "\t" + r.end);
    //  }
        return regions;
    }
    private void merge(){
        boolean merged = true;
        ArrayList<Region> buf = new ArrayList<>(this.regions);
        Collections.sort(buf);
        if(buf.size() == 0){
            System.err.println("No regions were loaded for analysis. Maybe gene names are bad or the file is broken.");
            System.exit(-1);
        }
        while(merged){
            merged = false;
            ArrayList<Region> pool = new ArrayList<>();
            pool.add(buf.get(0));
            // System.err.println("--------------");
            for(int i = 1; i<buf.size(); i++){
                Region last = pool.get(pool.size() -1 );
                if(buf.get(i).hasIntersection(last)){
                // System.out.println("has intersect with last record");
                // System.err.println("merged");
                   last.merge(buf.get(i));
                   merged = true;
                }else {
                   pool.add(buf.get(i));
                }
            }
            buf = pool;
        }
        this.regions = buf;
    }
    public static void main(String[] argv){
        try {
            ArrayList<String> genes = new ArrayList<String>();
            genes.add("ANKRD11");
            genes.add("NFKB1");
            genes.add("PADI4");
            File file = new File("resources/S04380110_Padded.bed");
            SureSelect ss = new SureSelect(file, genes, true);
            for(Region r: ss.getRegions()){
                System.out.println(r.gene + "\t" + r.chr + "\t" + r.start + "\t" + r.end);
            }
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
