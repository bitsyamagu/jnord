package jnord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;


public class Main {
    SureSelect ss;
    List<String> geneNames;
    ArrayList<Gene> genes = new ArrayList<>();
    ArrayList<String> samples = new ArrayList<>();
    HashSet<String> targetSamples = new HashSet<>();
    HashMap<String, ArrayList<RefGene>> refgeneMap = null;

    public Main(List<String> geneNames){
        this.geneNames = geneNames;
    }
    public void init(ArrayList<String> targetSamples) throws FileNotFoundException, IOException{
         File file = new File(AnalysisContext.captureBed);
         this.targetSamples.addAll(targetSamples);

         this.ss = new SureSelect(file, this.geneNames, true);
         // check and abort
         boolean hasError = false;
         for(String geneName: this.geneNames){
             if(this.ss.getRegions(geneName) == null){
                 System.err.println("Bad Gene name: " + geneName + "\n Couldn't find in this SureSelect");
                 hasError = true;
             }else if(this.ss.getRegions(geneName).size() == 0){
                 System.err.println("Bad Gene : " + geneName + "\n Insufficient gene information (no cds or no exon)");
                 hasError = true;
             }
         }
         if(hasError){
                 System.exit(0);
         }
         // System.err.println("start check===============");
         for(String geneName: this.geneNames){
             // System.err.println(geneName);
             Gene gene = new Gene(geneName, this.ss.getRegions(geneName));
             // System.err.println(geneName + "\t" + gene.getFirst().start + "\t" + gene.getLast().end);
             genes.add(gene);
             // gene could have multiple NM_* accession
         }
    }
    public void prepGCPercent(String path, String indexPath) throws FileNotFoundException{
        GCPercent gcp = new GCPercent(new File(path), new File(indexPath));
        for(Gene gene: this.genes){
            for(Region region: gene.getRegions()){
                region.setGCContent(gcp.getGCpercent100(region));
            }
        }
    }
    public void loadDepth(String bamDir) throws IOException{
        File dir = new File(bamDir);
        File[] bams = dir.listFiles(
            new java.io.FilenameFilter(){
                public boolean accept(File f, String filename){
                    return filename.endsWith(".bam");
                }
            }
        );
        HashMap<String, File> sampleBams = new HashMap<String, File>();
        for(File bam: bams){
            SamReader reader = SamReaderFactory.makeDefault().open(bam);
            SAMFileHeader h = reader.getFileHeader();
            List<SAMReadGroupRecord> rgrList = h.getReadGroups();
            for(SAMReadGroupRecord rgr: rgrList){
                String sampleName = rgr.getSample();
                samples.add(sampleName);
                System.err.println("loading depth of " + sampleName);
                sampleBams.put(sampleName, bam);
                for(Gene gene: this.genes){
                    for(Region region: gene.getRegions()){
                        IntervalList intervalList = new IntervalList(h);
                        intervalList.add(new Interval(region.getChromosome(), region.getStart(), region.getEnd()));
                        SamLocusIterator slit = new SamLocusIterator(reader, intervalList);
                        double[] buf = new double[region.getEnd() - region.getStart() + 1];
                        int i = 0;
                        // System.out.println("chr" + region.getChromosome() + " " + (region.getStart()));
                        while(slit.hasNext()){
                            SamLocusIterator.LocusInfo info = slit.next();
                            buf[i] = info.size();
                            i++;
                        }
                        slit.close();
                        region.addDepth(sampleName, buf);
                    }
                }
            }
            reader.close();
        }
    }
    public void calcMedian(){
        Percentile perc = new Percentile();
        for(Gene gene: this.genes){
            System.err.println("calcurating median for gene " + gene.getName());
            for(Region region: gene.getRegions()){
                ArrayList<double[]> depths = new ArrayList<>(region.getDepths());
                double[] median = new double[region.getLength()];
                double[] buf = new double[depths.size()];
                double[] sd = new double[region.getLength()];
                double[] sn = new double[region.getLength()];
                for(int i = 0; i<median.length; i++){
                    for(int j = 0; j<depths.size(); j++){
                        buf[j] = depths.get(j)[i];
                    }
                    Arrays.sort(buf);
                    median[i] = perc.evaluate(buf, 50.0f);
                    sd[i] = sd(buf);
                    sn[i] = median[i]/sd[i];
                }
                region.setMedian(median);
                region.setSD(sd);
                region.setSN(sn);
            }
        }
    }
    private double sd(double[] buf){
        double mean = 0.0f;
        for(double x: buf){
            mean += x;
        }
        mean = mean/buf.length;
        double v = 0.0f;
        for(double x: buf){
            v = v + (x-mean)*(x-mean);
        }
        return Math.sqrt(v/(double)buf.length);
    }
    public static List<String> readList(File f) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(f));
        String raw = null;
        ArrayList<String> list = new ArrayList<>();
        while(null != (raw = br.readLine())){
            if(raw.startsWith("#")){
                continue;
            }
            list.add(raw.trim());
        }
        br.close();
        return list;
    }
    public void normalize(){
        // count data point for autosome
        System.err.println("==== Normalization by InvariantSet ====");
        int count = 0;
        for(Gene gene: this.genes){
            for(Region r: gene.getRegions()){
                String chr = r.getChromosome();
                if(!(chr.startsWith("X") || chr.startsWith("Y") || chr.startsWith("M"))){
                    count += r.getLength();
                }
            }
        }
        double[] rawIntensities = new double[count]; // buffer for normalization
        double[] medians = new double[count];        // buffer for normalization
        for(String sample: this.samples){
            System.err.println("Normalizing for " + sample);
            count = 0;
            // collect values for mormalization on autosomes
            for(Gene gene: this.genes){
                for(Region r: gene.getRegions()){
                    String chr = r.getChromosome();
                    if(!(chr.startsWith("X") || chr.startsWith("Y") || chr.startsWith("M"))){
                        double[] depth = r.getDepth(sample);
                        double[] median = r.getMedian();
                        for(int i = 0; i<depth.length; i++){
                            rawIntensities[count] = depth[i];
                            medians[count] = median[i];
                            count++;
                        }
                    }
                }
            }
            // System.out.println(sample);
            InvariantSetNormalization normalizer = new InvariantSetNormalization();
            normalizer.init(rawIntensities, medians); // build spline function for normalization
            for(Gene gene: this.genes){
                for(Region r: gene.getRegions()){
                    double[] depth = r.getDepth(sample);
                    double[] median = r.getMedian();
                    double[] normalized = new double[depth.length];
                    double[] ratioNormalized = new double[depth.length];
                    for(int i = 0; i<depth.length; i++){
                        normalized[i] = normalizer.normalize(depth[i], median[i]);
                        ratioNormalized[i] = normalized[i] / median[i];
                    }
                    r.addNormalized(sample, normalized);
                    r.addRatioNormalized(sample, ratioNormalized);
                    /*
                    if(sample.equals("Sample_17063") && r.getStart() > 41545000 && r.getEnd() < 41547000){
                        for(int i = 0; i<depth.length; i++){
                            normalizer.setDebug(true);
                            normalized[i] = normalizer.normalize(depth[i], median[i]);
                            normalizer.setDebug(false);
                            ratioNormalized[i] = normalized[i] / median[i];
                            System.out.println("# raw_depth: " + depth[i] + "\tmedian: " 
                              + median[i] + "\tnorm: " + normalized[i]  + "\tratio: " + ratioNormalized[i]);
                        }
                    }*/
                }
            }
            if(AnalysisContext.plotNormalization){
                InvariantSetNormalizationPlot plotter = new InvariantSetNormalizationPlot(sample);
                try {
                    plotter.plot(this.genes, normalizer.getCurve());
                }catch(Exception e){
                    System.err.println("Cannot write plot for sample "+ sample);
                    throw new RuntimeException(e);
                }
            }
        }

        System.err.println("==== Normalization finished ====");
        // normalizer.init(prog.getAutosomeIntensities(), prog.getAutosomeMedians());
    }

    public void calcQC(){
        Percentile perc = new Percentile();
        for(Gene gene: this.genes){
            System.err.println("calcurating QC for gene " + gene.getName());
            for(Region region: gene.getRegions()){
                ArrayList<double[]> depths = new ArrayList<>(region.getNormalized());
                double[] buf = new double[depths.size()];
                double[] sn = new double[region.getLength()];
                for(int i = 0; i<sn.length; i++){
                    for(int j = 0; j<depths.size(); j++){
                        buf[j] = depths.get(j)[i];
                    }
                    Arrays.sort(buf);
                    sn[i] = perc.evaluate(buf, 50.0f)/sd(buf);
                }
                region.setNormalizedSN(sn);
            }
        }
    }
    public void callCNV(){
        Percentile perc = new Percentile();
        HashMap<String, ArrayList<CNV>> result = new HashMap<>();
        for(String sample: this.samples){
            result.put(sample, new ArrayList<CNV>());
        }
        for(Gene gene: this.genes){
            CNVCaller caller = new CNVCaller();
            caller.run(gene, samples);
            HashMap<String, ArrayList<CNV>> map = caller.getCNVs();
            // only target sampels are added
            for(String sample: this.targetSamples){
                result.get(sample).addAll(map.get(sample));
            }
        }
        // output target samples
        for(String sample: this.targetSamples){
            ArrayList<CNV> sampleCNV = result.get(sample);
            for(CNV cnv: sampleCNV){
                System.out.println(cnv.getChromosome() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + cnv.getType() + "\t" + sample+ "\t" + cnv.getGene().getName());
            }
            // plot
            HashMap<String, Gene> map = new HashMap<>();
            for(Gene gene: this.genes){
                map.put(gene.getName(), gene);
            }
            for(String geneName: this.geneNames){
                Gene gene = map.get(geneName);
                ArrayList<RefGene> targetGene = this.refgeneMap.get(gene.getName());
                // System.out.println(gene.getName() + " has " + targetGene.size() + " isoforms");

                GenePlot plot = new GenePlot(sample, targetGene, 1900, 600);
                for(Region region: gene.getRegions()){
                    plot.addData(region);
                }
                for(CNV cnv: sampleCNV){
                    if(cnv.getGene().getName().equals(gene.getName())){
                        plot.addCNV(cnv);
                    }
                }

                plot.paint();
                try {
                    plot.write(sample + "_" + gene.getName() + ".png", "png");
                }catch(Exception e){
                    e.printStackTrace();
                }
            }
        }
    }
    private void initPlot(String refgenePath){
        try {
            // ArrayList<RefGene> refgenes = RefGene.load(new java.util.zip.GZIPInputStream(new java.io.FileInputStream("./resources/UCSC-20181203/refGene.txt.gz")));
            ArrayList<RefGene> refgenes = RefGene.load(new java.util.zip.GZIPInputStream(new java.io.FileInputStream(refgenePath)));
            HashMap<String, ArrayList<RefGene>> map = new HashMap<>();
            for(RefGene refgene: refgenes){
                if(!map.containsKey(refgene.getName())){
                    map.put(refgene.getName(), new ArrayList<RefGene>());
                }
                map.get(refgene.getName()).add(refgene);
            }
            this.refgeneMap = map;
        }catch(Exception e){
            e.printStackTrace();
        }
    }
    public void printTable() throws IOException{
        TableWriter writer = new TableWriter(this.geneNames, this.genes, this.samples);
        writer.write();
    }
    public void drawAllGenesPlot() throws IOException, FileNotFoundException{
        ArrayList<String> plotSamples = new ArrayList<>(this.targetSamples);
        if(AnalysisContext.plotAllGenePlotForAllSamples){
            plotSamples = this.samples;
        }
        for(String sample: plotSamples){
            System.err.println("Drawing all gene plots for " + sample);
            AllGenePlot allGenePlot = new AllGenePlot(sample, true);
            allGenePlot.plot(this.genes);
        }
    }
      
    public static void main(String[] argv){
        String geneList = argv[0];
        String bamDir = argv[1];
        String refgenePath = AnalysisContext.refgene;
        ArrayList<String> targetSamples = new ArrayList<>();
        try {
            for(int i = 2; i<argv.length; i++){
                targetSamples.add(argv[i]);
            }
            Main prog = new Main(readList(new File(geneList)));
            prog.init(targetSamples);
            prog.initPlot(refgenePath);
            prog.loadDepth(bamDir);
            prog.calcMedian();
            // prog.prepGCPercent(argv[2], argv[3]);
            prog.normalize();
            prog.calcQC();
            prog.callCNV();
            prog.printTable();
            prog.drawAllGenesPlot();

        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
