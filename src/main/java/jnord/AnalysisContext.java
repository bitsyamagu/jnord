package jnord;

public class AnalysisContext {
    // SureSelect Capture information: padded file includes outside bases of target exons
    // public static String captureBed = "/data/jNord/resources/S04380110_Padded.bed";
    // public static String refgene = "/data/jNord/resources/UCSC-20181203/refGene.txt.gz";
    public static String captureBed = "/usr/local/bio/src/jnord/resources/S04380110_Padded.bed";
    public static String refgene = "/usr/local/bio/src/jnord/resources/refGene.txt.gz";

    // reference sequeces are required for GC percent normalization(not in use 201812)
    public static String reference = "/usr/local/bio/db/riker/bundle_b37_novoalign/human_g1k_v37_fix.fasta";
    public static String referenceIndex = "/usr/local/bio/db/riker/bundle_b37_novoalign/human_g1k_v37_fix.fasta.fai";

    // Image output optoins
    public static boolean plotAllGenePlotForAllSamples = false;
    public static boolean plotNormalization = true;
    
    // CNV Caller configuration
    public static double[] cnvCallerLossDetectionThresholds = {0.7, 0.8};
    public static double[] cnvCallerGainDetectionThresholds = {1.3, 1.2};
    public static double[] cnvCallerNoneDetectionThresholds = {0.2, 0.3};
    public static double[] cnvCallerManyDetectionThresholds = {1.8, 1.7};
    // private Detector lossDetector = new LessDetector(0.7, 0.8, minCoverage);
    // private Detector gainDetector = new MoreDetector(1.3, 1.2, minCoverage);
    // private Detector lossNoneDetector = new LessDetector(0.2, 0.3, minCoverage);
    // private Detector gainManyDetector = new MoreDetector(1.8, 1.7, minCoverage);
    public static int cnvCallerMinCoverage = 30;
    public static int cnvCallerWindowSize = 50;

    public static boolean randomColorScheme = false;
    public static boolean legend = true;

    public static double SN_threshold = 3.0f;

    public static boolean debug = false;
}
