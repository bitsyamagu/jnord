package jnord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;


/*
 * multi-threaded version of Main
 */
public class Main2 {
    SureSelect ss;
    List<String> geneNames;
    ArrayList<Gene> genes = new ArrayList<>();
    ArrayList<String> samples = new ArrayList<>();
    HashSet<String> targetSamples = new HashSet<>();
    HashMap<String, ArrayList<RefGene>> refgeneMap = null;
    ThreadPoolExecutor executor = null;
    int threadPoolSize = 1;

    public Main2(List<String> geneNames){
        this.geneNames = geneNames;
    }
    private void prepareThreadPool(){
         this.executor = (ThreadPoolExecutor)Executors.newFixedThreadPool(threadPoolSize);
    }
    public void setSureSelect(String path){
        AnalysisContext.captureBed = path;
    }
    public void init(ArrayList<String> targetSamples, int threadpoolSize) throws FileNotFoundException, IOException{
         this.threadPoolSize = threadpoolSize;
         File file = new File(AnalysisContext.captureBed);
         this.targetSamples.addAll(targetSamples);
         boolean hasError = false;

         if(this.geneNames.size() < 1){
             System.err.println("No genes in the file");
             System.exit(0);
         }
         this.ss = new SureSelect(file, this.geneNames, true);
         // check and abort
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
    public ArrayList<String> loadSamples(String bamDir) throws IOException {
        File dir = new File(bamDir);
        if(!dir.exists()){
            System.err.println("BAM directory " + dir.getName() + " doesn't exist");
            System.exit(1);
        }
        File[] bams = dir.listFiles(
            new java.io.FilenameFilter(){
                public boolean accept(File f, String filename){
                    return filename.endsWith(".bam");
                }
            }
        );
        HashSet<String> samples = new HashSet<String>();
        // Executor
        for(File bam: bams){
            SamReader reader = SamReaderFactory.makeDefault().open(bam);
            SAMFileHeader h = reader.getFileHeader();
            List<SAMReadGroupRecord> rgrList = h.getReadGroups();
            for(SAMReadGroupRecord rgr: rgrList){
                String sampleName = rgr.getSample();
                samples.add(sampleName);
                System.err.println("Sample " + sampleName + " will be loaded for analysis");
            }
            reader.close();
        }
        return new ArrayList<String>(samples);
    }
    public void loadDepth(String bamDir) throws IOException{
        File dir = new File(bamDir);
        if(!dir.exists()){
            System.err.println("BAM directory " + dir.getName() + " doesn't exist");
            System.exit(1);
        }
        File[] bams = dir.listFiles(
            new java.io.FilenameFilter(){
                public boolean accept(File f, String filename){
                    return filename.endsWith(".bam");
                }
            }
        );
        HashMap<String, File> sampleBams = new HashMap<String, File>();
        // Executor
        prepareThreadPool();
        for(File bam: bams){
            this.executor.execute( () -> {
              try {
                SamReader reader = SamReaderFactory.makeDefault().open(bam);
                SAMFileHeader h = reader.getFileHeader();
                List<SAMReadGroupRecord> rgrList = h.getReadGroups();
                for(SAMReadGroupRecord rgr: rgrList){
                    String sampleName = rgr.getSample();
                    synchronized(this.samples){
                        this.samples.add(sampleName);
                        System.err.println("loading depth of " + sampleName);
                    }
                    synchronized(sampleBams){
                        sampleBams.put(sampleName, bam);
                    }
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
                    // checking run information
                    SAMRecordIterator rit = reader.iterator();
                    SAMRecord rec = rit.next();
                    System.err.println("Readname for " + sampleName + ": " + rec.getReadName());
                }
                reader.close();
              }catch(IOException e){
                  System.err.println("I/O Error at " + bam.getPath() + " "+ e.getMessage());
              }
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
        }
    }
    public void calcMedian(){
        // Executor
        prepareThreadPool();
        for(Gene gene: this.genes){
            this.executor.execute( () -> {
                Percentile perc = new Percentile();
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
                        if(Double.isInfinite(sn[i])){
                            sn[i] = 0.0;
                        }
                    }
                    region.setMedian(median);
                    region.setSD(sd);
                    region.setSN(sn);
                }
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
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
        if(!f.exists()){
            System.err.println("The file "+ f.getName() + " does not exist");
            System.exit(-1);
        }
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
        int length = 0;
        for(Gene gene: this.genes){
            for(Region r: gene.getRegions()){
                String chr = r.getChromosome();
                if(!(chr.startsWith("X") || chr.startsWith("Y") || chr.startsWith("M"))){
                    length += r.getLength();
                }
            }
        }
        if(length == 0){
            throw new RuntimeException("Please include autosomes in this analysis for normalization");
        }
        prepareThreadPool();
        final int len = length;
        for(String sample: this.samples){
            this.executor.execute( () -> {
                System.err.println("Normalizing for " + sample);
                double[] rawIntensities = new double[len]; // buffer for normalization
                double[] medians = new double[len];        // buffer for normalization
                int count = 0;
                // collect values for mormalization on autosomes
                for(Gene gene: this.genes){
                    for(Region r: gene.getRegions()){
                        String chr = r.getChromosome();
                        if(!(chr.startsWith("X") || chr.startsWith("Y") || chr.startsWith("M") || chr.startsWith("GL"))){
                            double[] depth = r.getDepth(sample);
                            double[] median = r.getMedian();
                            if(depth == null){
                                System.err.println("bad data for " + sample + " at " + chr + " for " + gene.getName());
                                System.exit(-1);
                            }
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
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
        }

        System.err.println("==== Normalization finished ====");
        // normalizer.init(prog.getAutosomeIntensities(), prog.getAutosomeMedians());
    }

    public void calcQC(){
        prepareThreadPool();
        for(Gene gene: this.genes){
            this.executor.execute( () -> {
                Percentile perc = new Percentile();
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
                        sn[i] = perc.evaluate(buf, 50.0f)/sd(buf); // (median of normalized depth)/SD
                    }
                    region.setNormalizedSN(sn);
                }
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
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
                if(map.get(sample) == null){
                    System.err.println("ERROR: The data for " + sample + " is not ready yet.");
                    System.exit(-1);
                }
                result.get(sample).addAll(map.get(sample));
            }
        }
        prepareThreadPool();
        // output target samples
        for(String sample: this.targetSamples){
            this.executor.execute( () -> {
                ArrayList<CNV> sampleCNV = result.get(sample);
                for(CNV cnv: sampleCNV){
                    System.out.println(cnv.getChromosome() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + cnv.getType() + "\t" + sample+ "\t" + cnv.getGene().getName());
                }
                // plot
                for(Gene gene: this.genes){
                    if(!this.refgeneMap.containsKey(gene.getName())){
                        continue;
                    }
                    ArrayList<RefGene> targetGene = this.refgeneMap.get(gene.getName());
                    // System.out.println(gene.getName() + " has " + targetGene.size() + " isoforms");

                 // if(targetGene == null){
                 //     System.err.println("The gene " + gene.getName() + " cannot find in RefGene");
                 // }
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
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
        }
    }
    private void checkGenes(){
        boolean hasError = false;
        for(Gene gene: this.genes){
            if(!this.refgeneMap.containsKey(gene.getName())){
                System.err.println("The gene " + gene.getName() + " does not exist in RefGene");
                hasError = true;
            }
        }
        if(hasError){
            // System.exit(1);
            System.err.println("These genes cannot be plotted");
        }
    }
    private void initPlot(String refgenePath){
        try {
            // ArrayList<RefGene> refgenes = RefGene.load(new java.util.zip.GZIPInputStream(new java.io.FileInputStream("./resources/UCSC-20181203/refGene.txt.gz")));
            ArrayList<RefGene> refgenes = null;
            if(refgenePath.endsWith(".gz")){
                refgenes = RefGene.load(new java.util.zip.GZIPInputStream(new java.io.FileInputStream(refgenePath)));
            }else {
                refgenes = RefGene.load(new java.io.FileInputStream(refgenePath));
            }
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
        java.util.Collections.sort(this.samples);
        TableWriter writer = new TableWriter(this.geneNames, this.genes, this.samples);
        writer.write();
    }
    public void drawAllGenesPlot() throws IOException, FileNotFoundException{
        java.util.Collections.sort(this.genes);
        ArrayList<String> plotSamples = new ArrayList<>(this.targetSamples);
        if(AnalysisContext.plotAllGenePlotForAllSamples){
            plotSamples = this.samples;
        }
        prepareThreadPool();
        for(String sample: plotSamples){
            this.executor.execute( () -> {
                try {
                    System.err.println("Drawing all gene plots for " + sample);
                    AllGenePlot allGenePlot = new AllGenePlot(sample, AnalysisContext.legend);
                    allGenePlot.plot(this.genes);
                }catch(IOException e){
                    System.err.println(e.getMessage());
                    e.printStackTrace();
                }
            });
        }
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE,TimeUnit.NANOSECONDS);
        }catch(InterruptedException ie){
            System.err.println(ie.getMessage());
            System.exit(-1);
        }
    }
      
    public static void main(String[] argv){
        String geneList = null;
        String bamDir = null;
        int threadpoolSize = 4;
        String sureselect = null;
        ArrayList<String> targetSamples = new ArrayList<>();
        boolean allSamples = false;
        for(int i = 0; i<argv.length-1; i++){
            if(argv[i].equals("--threads")){
                threadpoolSize = Integer.parseInt(argv[i+1]);
            }else if(argv[i].equals("--plot-width")){
                AnalysisContext.auto_scale = false;
                AnalysisContext.plot_width = Integer.parseInt(argv[i+1]);
            }else if(argv[i].equals("--bamdir")){
                bamDir = argv[i+1];
            }else if(argv[i].equals("--samples")){
                if(argv[i+1].equals("all")){
                    allSamples = true;
                }
                for(String sample: argv[i+1].split(",")){
                    targetSamples.add(sample);
                }
            }else if(argv[i].equals("--genes")){
                geneList = argv[i+1];
            }else if(argv[i].equals("--cnvCallerMinCoverage")){
                AnalysisContext.cnvCallerMinCoverage = Integer.parseInt(argv[i+1]);
            }else if(argv[i].equals("--sn_threshold")){
                AnalysisContext.SN_threshold = Double.parseDouble(argv[i+1]);
            }else if(argv[i].equals("--sureselect")){
                sureselect = argv[i+1];
            }else if(argv[i].equals("--refgene")){
                AnalysisContext.refgene = argv[i+1];
            }else if(argv[i].equals("--randomcolor")){
                AnalysisContext.randomColorScheme = Boolean.parseBoolean(argv[i+1]);
            }else if(argv[i].equals("--legend")){
                if(argv[i+1].equals("true")){
                    AnalysisContext.legend = true;
                }else if(argv[i+1].equals("false")) {
                    AnalysisContext.legend = false;
                }else {
                    System.err.println("Unknown option for --legend");
                    System.exit(-1);
                }
            }else if(argv[i].equals("--debug")){
                AnalysisContext.debug = Boolean.parseBoolean(argv[i+1]);
            }else if(argv[i].equals("--windowfilter")){
                // --windowfilter 180,200
                String[] thresholds = argv[i+1].split(",");
                AnalysisContext.windowFilter_passThreshold = Integer.parseInt(thresholds[0]);
                AnalysisContext.windowFilter_lengthThreshold = Integer.parseInt(thresholds[1]);
            }else if(argv[i].startsWith("-")){
                System.err.println("Unrecognized option " + argv[i]);
                System.exit(-1);
            }
        }
        String refgenePath = AnalysisContext.refgene;
        boolean hasError = false;
        if(geneList == null){
            System.err.println("An gene list file is mandatory for analysis: use --genes [file] option");
            hasError = true;
        }
        if(!allSamples && targetSamples.size() == 0){
            System.err.println("Use --samples [sample1,sample2,sample3...] option to identify target samples.");
            hasError = true;
        }
        if(bamDir == null){
            System.err.println("Use --bamdir [dir] option for the directory containing target bam files");
            hasError = true;
        }

        if(hasError){
            System.exit(-1);
        }
        try {
         // for(int i = 2; i<argv.length; i++){
         //     targetSamples.add(argv[i]);
         // }
            Main2 prog = new Main2(readList(new File(geneList)));
            if(sureselect != null){
                prog.setSureSelect(sureselect);
            }
            if(allSamples){
                targetSamples = prog.loadSamples(bamDir);
                java.util.Collections.sort(targetSamples);
            }
            prog.init(targetSamples, threadpoolSize);
            prog.initPlot(refgenePath);
            prog.checkGenes();
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
