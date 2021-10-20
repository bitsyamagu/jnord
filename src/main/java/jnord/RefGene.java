package jnord;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.zip.GZIPInputStream;

public class RefGene extends Gene{
    String accession;
    String chr;
    String strand;
    int transcriptionStart;
    int transcriptionEnd;
    int cdsStart;
    int cdsEnd;
    Exon[] exons = null;
    String score;
    String cdsStartStat;
    String cdsEndStat;
    String exonFrames;

    public RefGene(String name){
        super(name);
    }
    public RefGene(String acc, String chr, String strand, int txStart, int txEnd, int cdsStart, int cdsEnd, Exon[] exons, String score, String name, String cdsStartStat, String cdsEndStat, String exonFrames){
        super(name);
        this.accession = acc;
        this.chr = chr;
        this.strand = strand;
        this.transcriptionStart = txStart;
        this.transcriptionEnd = txEnd;
        this.cdsStart = cdsStart;
        this.cdsEnd = cdsEnd;
        this.exons = exons;
        this.score = score;
        this.name = name;
        this.cdsStartStat = cdsStartStat;
        this.cdsEndStat = cdsEndStat;
        this.exonFrames = exonFrames;
    }
    public Exon getFirstExon(){
        return this.exons[0];
    }
    public Exon getLastExon(){
        return this.exons[exons.length-1];
    }
    public static ArrayList<RefGene> load(InputStream is) throws IOException{
        ArrayList<RefGene> genes = new ArrayList<>();
        BufferedReader br = new BufferedReader(new InputStreamReader(is));
        String raw = null;
        while(null != (raw = br.readLine())){
            String[] line = raw.split("\t", 16);
            if(line.length < 16){
                System.err.println("skipping: " + raw);
                continue;
            }
            String acc = line[1];
            String chr = line[2];
            if(chr.startsWith("chr")){
                chr = chr.substring(3);
            }
            String strand = line[3];
            int txStart = Integer.parseInt(line[4]);
            int txEnd = Integer.parseInt(line[5]);
            int cdsStart = Integer.parseInt(line[6]);
            int cdsEnd   = Integer.parseInt(line[7]);
            int exonCount = Integer.parseInt(line[8]);
            Exon[] exons = new Exon[exonCount];
            String[] exon_starts = line[9].split(",");
            String[] exon_ends = line[10].split(",");
            for(int i = 0; i<exonCount; i++){
                exons[i] = new Exon(chr, Integer.parseInt(exon_starts[i]), Integer.parseInt(exon_ends[i]));
            }
            String score = line[11];
            String name = line[12];
            String cdsStartStat = line[13];
            String cdsEndStat = line[14];
            String exonFrame = line[15];
            RefGene rg = new RefGene(acc, chr, strand, txStart, txEnd, cdsStart, cdsEnd, exons, score, name, cdsStartStat, cdsEndStat, exonFrame);
            genes.add(rg);
        }
        br.close();
        return genes;
    }
    public static ArrayList<Exon> mergeExons(ArrayList<RefGene> genes){
        ArrayList<Exon> exons = new ArrayList<>();
        for(RefGene g: genes){
            for(Exon ex: g.getExons()){
                exons.add(ex);
            }
        }
        java.util.Collections.sort(exons);
        boolean flag = true;
        while(flag){
            flag = false;
            INNER:
            for(int i = 1; i<exons.size(); i++){
                Exon prev = exons.get(i-1);
                Exon cur = exons.get(i);
                if(prev.hasIntersection(cur)){
                    exons.remove(cur);
                    prev.merge(cur);
                    flag = true;
                    break INNER;
                }
            }
        }
   //   System.err.println("--merged--");
   //   for(Exon exon: exons){
   //       System.err.println(exon.getStart() + " "+ exon.getEnd());
   //   }
   //   System.err.println("----");
        return exons;
    }
    public Exon[] getExons(){
        return this.exons;
    }

    public static void main(String[] argv){
        String path = "resources/UCSC-20181203/refGene.txt.gz";
        try {
            ArrayList<RefGene> genes = RefGene.load(new GZIPInputStream(new FileInputStream(path)));
            HashMap<String, ArrayList<RefGene>> map = new HashMap<>();
            for(RefGene gene: genes){
                if(!map.containsKey(gene.getName())){
                    map.put(gene.getName(), new ArrayList<RefGene>());
                }
                map.get(gene.getName()).add(gene);
            }
    
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
// 0    1              2       3   4           5           6           7           8   9                                                                                           10                                                                                                                                                  11  12      13      14       15
// 158  NR_045839      chr16   -   89353065    89556969    89556969    89556969    10  89353065,89354935,89357032,89357420,89358088,89371613,89379724,89383340,89484691,89556652,  89353521,89355078,89357236,89357591,89358185,89371752,89379997,89383486,89484776,89556969,                                                          0   ANKRD11 unk     unk     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
// 158  NM_001256183   chr16   -   89334028    89556969    89334885    89383427    13  89334028,89337224,89341221,89341500,89345479,89352446,89354935,89357032,89357420,89371613,89383340,89484691,89556652,   89335071,89337317,89341365,89341599,89352057,89352594,89355078,89357236,89357591,89371752,89383483,89484776,89556969,   0   ANKRD11 cmpl    cmpl    0,0,0,0,1,0,1,1,1,0,0,-1,-1,
// 158  NM_001256182   chr16   -   89334028    89556969    89334885    89383427    14  89334028,89337224,89341221,89341500,89345479,89352446,89354935,89357032,89357420,89371613,89383340,89484691,89497663,89556652,  89335071,89337317,89341365,89341599,89352057,89352594,89355078,89357236,89357591,89371752,89383486,89484776,89497734,89556969,  0   ANKRD11 cmpl    cmpl    0,0,0,0,1,0,1,1,1,0,0,-1,-1,-1,
// 158  NM_013275      chr16   -   89334028    89556969    89334885    89383427    13  89334028,89337224,89341221,89341500,89345479,89352446,89354935,89357032,89357420,89371613,89383340,89484691,89556652,   89335071,89337317,89341365,89341599,89352057,89352594,89355078,89357236,89357591,89371752,89383486,89484776,89556969,   0   ANKRD11 cmpl    cmpl    0,0,0,0,1,0,1,1,1,0,0,-1,-1,
//
