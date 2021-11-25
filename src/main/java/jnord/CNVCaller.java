package jnord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Collections;

public class CNVCaller {
    final static int CORE = 2;
    final static int EXTEND = 1;
    final static int NONE = 0;
    final static int CNV_MERGE_THRESHOLD = 5000;

    public class LessDetector implements Detector {
        double core;
        double extend;
        double minMedian;
        public LessDetector(double core_, double extend_, double minMedian_){
            core = core_;
            extend = extend_;
            this.minMedian = minMedian_;
        }
        public int check(double score, double median){
            if(score < core && median >= minMedian){
                return CORE;
            }else if(score < extend && median >= minMedian){
                return EXTEND;
            }
            return NONE;
        }
    };
    public class MoreDetector implements Detector {
        double core;
        double extend;
        double minMedian;
        public MoreDetector(double core_, double extend_, double minMedian_){
            core = core_;
            extend = extend_;
            minMedian = minMedian_;
        }
        public int check(double score, double median){
            if(score > core && median >= minMedian){
                return CORE;
            }else if(score > extend && median >= minMedian){
                return EXTEND;
            }
            return NONE;
        }
    };
    // default detectors
    double minCoverage = AnalysisContext.cnvCallerMinCoverage; // default 50
    private Detector lossDetector = new LessDetector(
        AnalysisContext.cnvCallerLossDetectionThresholds[0],
        AnalysisContext.cnvCallerLossDetectionThresholds[1], minCoverage);
    private Detector gainDetector = new MoreDetector(
        AnalysisContext.cnvCallerGainDetectionThresholds[0],
        AnalysisContext.cnvCallerGainDetectionThresholds[1], minCoverage);
    private Detector lossNoneDetector = new LessDetector(
        AnalysisContext.cnvCallerNoneDetectionThresholds[0],
        AnalysisContext.cnvCallerNoneDetectionThresholds[1], minCoverage);
    private Detector gainManyDetector = new MoreDetector(
        AnalysisContext.cnvCallerManyDetectionThresholds[0],
        AnalysisContext.cnvCallerManyDetectionThresholds[1], minCoverage);

    private int minHitBaseCount = 37;
    HashMap<String, ArrayList<CNV>> map = new HashMap<String, ArrayList<CNV>>();

    public CNVCaller(){
    }
    public HashMap<String, ArrayList<CNV>> getCNVs(){
        return this.map;
    }
    public void run(Gene gene, ArrayList<String> samples){
        int windowSize = AnalysisContext.cnvCallerWindowSize;
        System.err.println("=== Running CNV Caller for "+ gene.getName() +" ===");
        for(String sample: samples){

            ArrayList<CNV> loss = new ArrayList<>();
            ArrayList<CNV> gain = new ArrayList<>();
            ArrayList<CNV> many = new ArrayList<>();
            ArrayList<CNV> none = new ArrayList<>();
            // for loss
            for(Region region: gene.getRegions()){
                int start = region.getStart();
                int end = region.getEnd();
                if(AnalysisContext.debug){
                    System.err.println("checking " +sample + "\t" + start + "\t" + end);
                }

                double[] ratio = region.getRatioNormalized(sample);
         //     if(sample.equals("Sample_17063") && gene.getName().equals("ARID1B") && region.getStart() > 157200000 && region.getStart() < 157400000){
         //         for(int i = 0; i<ratio.length; i++){
         //             System.out.println("Sample_17063.ARID1B: " + ratio[i] + " at "+ (region.getStart() + i));
         //         }
         //     }
                double[] median = region.getMedian();
                if(ratio == null){
                    System.err.println("Error: ratio for "+  sample + " is null");
                    System.exit(-1);
                }
                int[] lossFlags = new int[ratio.length];
                int[] gainFlags = new int[ratio.length];
                int[] noneFlags = new int[ratio.length];
                int[] manyFlags = new int[ratio.length];
                
                for(int i = 0; i<ratio.length - windowSize; i++){
                    lossFlags[i] = lossDetector.check(ratio[i], median[i]);
                    gainFlags[i] = gainDetector.check(ratio[i], median[i]);
                    noneFlags[i] = lossNoneDetector.check(ratio[i], median[i]);
                    manyFlags[i] = gainManyDetector.check(ratio[i], median[i]);
                }
                ArrayList<Segment> lossSeeds = buildSeed(lossFlags);
                ArrayList<Segment> gainSeeds = buildSeed(gainFlags);
                ArrayList<Segment> noneSeeds = buildSeed(noneFlags);
                ArrayList<Segment> manySeeds = buildSeed(manyFlags);
                for(Segment s: lossSeeds){
                    loss.add(createCNV(s, region, "LOSS", gene));
                    if(AnalysisContext.debug){
                        System.err.println(sample + "\t" + region.getChromosome() + "\t " + region.getStart()  + "\t" + region.getEnd() + "\t" + "LOSS");
                    }
                }
                for(Segment s: gainSeeds){
                    gain.add(createCNV(s, region, "GAIN", gene));
                    if(AnalysisContext.debug){
                        System.err.println(sample + "\t" + region.getChromosome() + "\t " + region.getStart()  + "\t" + region.getEnd() + "\t" + "GAIN");
                    }
                }
                for(Segment s: manySeeds){
                    many.add(createCNV(s, region, "MANY", gene));
                    if(AnalysisContext.debug){
                        System.err.println(sample + "\t" + region.getChromosome() + "\t " + region.getStart()  + "\t" + region.getEnd() + "\t" + "MANY");
                    }
                }
                for(Segment s: noneSeeds){
                    none.add(createCNV(s, region, "NONE", gene));
                    if(AnalysisContext.debug){
                        System.err.println(sample + "\t" + region.getChromosome() + "\t " + region.getStart()  + "\t" + region.getEnd() + "\t" + "NONE");
                    }
                }
           //   if(manySeeds.size() > 0 || gainSeeds.size() > 0 || lossSeeds.size() > 0 || noneSeeds.size() > 0){
           //       System.out.println("Loss segments: " + lossSeeds.size());
           //       System.out.println("Gain segments: " + gainSeeds.size());
           //       System.out.println("None segments: " + noneSeeds.size());
           //       System.out.println("Many segments: " + manySeeds.size());
           //   }
            }
            Collections.sort(loss);
            Collections.sort(gain);
            Collections.sort(many);
            Collections.sort(none);
            // merge seeds
            loss = mergeCNV(loss);
            gain = mergeCNV(gain);
            many = mergeCNV(many);
            none = mergeCNV(none);

            loss = filter(sample, loss, lossDetector);
            gain = filter(sample, gain, gainDetector);
            none = filter(sample, none, lossNoneDetector);
            many = filter(sample, many, gainManyDetector);

            loss = filterLossOrGain(loss, none);
            gain = filterLossOrGain(gain, many);
            ArrayList<CNV> cnvs = new ArrayList<CNV>();
            cnvs.addAll(loss);
            cnvs.addAll(gain);
            cnvs.addAll(many);
            cnvs.addAll(none);
            Collections.sort(cnvs);
            this.map.put(sample, cnvs);
            // System.out.println(sample);
        }
    }
    private ArrayList<CNV> filter(String sample, ArrayList<CNV> cnvs, Detector detector){
        // filter by Median and Criteria
        MedianFilter medianFilter = new MedianFilter();
        WindowFilter windowFilter = new WindowFilter(AnalysisContext.windowFilter_passThreshold , AnalysisContext.windowFilter_lengthThreshold);
        ArrayList<CNV> result = new ArrayList<CNV>();
        for(CNV cnv: cnvs){
            double[] buf = cnv.getRatioNormalized(sample, cnv.getStart(), cnv.getEnd());
            double[] median = cnv.getMedian(cnv.getStart(), cnv.getEnd());
            if(medianFilter.accept(buf) && windowFilter.accept(buf, median, detector)){
                result.add(cnv);
                // System.err.print("PASS:\t"+ sample + "\t" + cnv.getChromosome() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t");
                // System.err.print("Median filter: " + medianFilter.accept(buf) + "\t" );
                // System.err.println("Window filter: " + windowFilter.accept(buf, median, detector) );
            } else if(AnalysisContext.debug) {
                System.err.print("FILT" + "\t" + detector.getClass().getName() + "\t" +sample + "\t" + cnv.getChromosome() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t");
                System.err.print("Median filter: " + medianFilter.accept(buf) + "\t" );
                System.err.println("Window filter: " + windowFilter.accept(buf, median, detector) );
            }
        }
        return result;
    }
    // filter loss and gain by none and many
    private ArrayList<CNV> filterLossOrGain(ArrayList<CNV> cnv1, ArrayList<CNV> cnv2){
        ArrayList<CNV> result = new ArrayList<CNV>();
        LOOP:
        for(CNV c1 : cnv1){
            boolean hit = false;
            for(CNV c2: cnv2){
                if(c2.hasIntersection(c1)){
                    hit = true;
                    continue LOOP;
                }
            }
            if(!hit){
                result.add(c1);
            }
            hit = false;
        }
        return result;
    }
    public ArrayList<CNV> mergeCNV(ArrayList<CNV> src){
        if(src.size() < 2){
            return src;
        }
        ArrayList<ArrayList<CNV>> list = new ArrayList<>();
        CNV prev = src.get(0);
        ArrayList<CNV> chain = new ArrayList<>();
        chain.add(prev);
        list.add(chain);
        for(int i = 1; i<src.size(); i++){
            CNV cur = src.get(i);
            if(canMerge(prev, cur)){
                chain.add(cur);
            }else {
                chain = new ArrayList<CNV>();
                chain.add(cur);
                list.add(chain);
            }
            prev = cur;
        }
        ArrayList<CNV> result = new ArrayList<>();
        for(ArrayList<CNV> c: list){
            CNV first = c.get(0);
            CNV last = c.get(c.size() - 1);
            result.add(new CNV(first.getChromosome(), first.getStart(), last.getEnd(), first.getType(), first.getGene()));
        }

        return result;
    }
    public CNV createCNV(Segment s, Region r, String cnvType, Gene gene){
        return new CNV(r.getChromosome(), r.getStart() + s.start, r.getStart() + s.end, cnvType, gene);
    }
    public boolean canMerge(CNV cnv1, CNV cnv2){
        if(!cnv1.getChromosome().equals(cnv2.getChromosome())){
            return false;
        }
   //   if(cnv1.getStart() > cnv2.getStart()){
   //       System.err.println("bad order");
   //       System.exit(1);
   //   }
        int spacing = cnv2.getStart() - cnv1.getEnd();
        boolean same = (cnv1.getStart() == cnv2.getStart() || cnv1.getEnd() == cnv2.getEnd());
        if(spacing < cnv1.getLength()*0.1 && spacing < cnv2.getLength()*0.1 && (spacing < CNV_MERGE_THRESHOLD || same)){
            return true;
        }
        return false;
    }
    private class Segment {
        int start;
        int end;
        public Segment(int s, int e) {
            this.start = s;
            this.end = e;
        }
    }
    private ArrayList<Segment> buildSeed(int[] flags){
        int windowSize = AnalysisContext.cnvCallerWindowSize;
        ArrayList<Integer> fillIndex = new ArrayList<>();
        for(int i = 0; i<flags.length - windowSize; i++){
            int cores = 0;
            int j = 0;
            for(j = i; j < i+windowSize; j++){
                if(flags[j] > 1){
                    cores++;
                }
            }
            if(cores > minHitBaseCount){
                // passed : cleared bases > 37 bases
                fillIndex.add(i);
            }
        }
        ArrayList<Segment> list = new ArrayList<Segment>();
        Segment lastSegment = null;
        for(Integer index: fillIndex){
            if(lastSegment == null){
                lastSegment = new Segment(index, index+windowSize); // initial state
            }else if(index < lastSegment.end){ // is overlapped
                lastSegment.end = index + windowSize; // update end
            }else {
                list.add(lastSegment);
            }
        }
        if((list.size() == 0 && lastSegment != null) || (list.size() > 0 && list.get(list.size()-1) != lastSegment)){
            list.add(lastSegment);
        }
        // fill and extend
        for(Segment s: list){
            // fill
            for(int i = s.start; i<s.end; i++){
                flags[i] = 2;
            }
            // extend for upstream
            for(int i = s.start; i>=s.start && flags[i]>0; i--){
                flags[i] = 2;
            }
            // extend for downstream
            for(int i = s.end; i<flags.length && flags[i] > 0; i++){
                flags[i] = 2;
            }
        }
        // fill if scores between segments are 1
        if(list.size() > 1){
            Segment cur = list.get(0);
            for(int i = 1; i<list.size(); i++){
                Segment next = list.get(i);
                for(int j = cur.end; j < next.start; j++){
                    if(flags[j] == 1){
                        flags[j] = 2;
                    }
                }
                cur = next;
            }
        }
        // rescan
        fillIndex = new ArrayList<>();
        for(int i = 0; i<flags.length - windowSize; i++){
            int cores = 0;
            int j = 0;
            for(j = i; j < i+windowSize; j++){
                if(flags[j] > 1){
                    cores++;
                }
            }
            if(cores > minHitBaseCount){
                // passed : cleared bases > 37 bases
                fillIndex.add(i);
            }
        }
        list = new ArrayList<Segment>();
        lastSegment = null;
        for(Integer index: fillIndex){
            if(lastSegment == null){
                lastSegment = new Segment(index, index+windowSize); // initial state
            }else if(index < lastSegment.end){ // is overlapped
                lastSegment.end = index + windowSize; // update end
            }else {
                list.add(lastSegment);
            }
        }
        if((list.size() == 0 && lastSegment != null) || (list.size() > 0 && list.get(list.size()-1) != lastSegment)){
            list.add(lastSegment);
        }
        //            System.out.println("//"+ cnv.getType() + " " + cnv.getChromosome() + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + cnv.getType() + "\t" + sample+ "\t" + cnv.getGene().getName());
        return list;
     }
}
