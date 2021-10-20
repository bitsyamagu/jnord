package jnord;

import java.awt.Dimension;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.Shape;
import java.awt.Insets;
import java.awt.Stroke;
import java.awt.BasicStroke;
import java.util.ArrayList;
import java.util.HashMap;
import java.awt.geom.Rectangle2D;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;

public class GenePlot {
    ArrayList<RefGene> gene = new ArrayList<>();
    int totalRegionLength = 0;
    Dimension size = new Dimension(800, 600);
    BufferedImage image = null;
    Insets insets = new Insets(20, 120, 80, 20); // top left bottom right
    HashMap<Integer, Integer> regionMap = new HashMap<>();
    ArrayList<Region> regions = new ArrayList<>();
    ArrayList<CNV> cnvs = new ArrayList<>();

    double REFGENE_TRACK_BASE_POS = 3.5;
    double GENOME_AXIS_POS = 3.65;
    int geneThickness = 10;
    int geneMargin = 3;
    int geneCount = 1;
    Color geneBorderColor = Color.black;
    Color geneBackgroundColor = Color.white;
    Color intronColor = Color.black;
    double logicalBottom = -2;
    double logicalTop = 4;
    RealMatrix convMatrix = null;
    int yoffset = 0;
    String sample;
    static final HashMap<String, Color> colorMap = new HashMap<>();
    static {
        colorMap.put("NONE", Color.pink);
        colorMap.put("LOSS", Color.red);
        colorMap.put("GAIN", Color.blue);
        colorMap.put("MANY", Color.cyan);
    };
    public GenePlot(String sample, ArrayList<RefGene> gene, int width, int height){
        this.sample = sample;
        this.gene = gene;
        this.setSize(width, height);
        // this.capturedRegion = RefGene.mergeExons(this.gene);
        // java.util.Collections.sort(capturedRegion);
    }
    public void setGene(ArrayList<RefGene> gene){
        this.gene = gene;
    }
    public void setSize(int width, int height){
        this.size = new Dimension(width, height);
    }
    public void addCNV(CNV cnv){
        this.cnvs.add(cnv);
    }
    public void addData(Region region){
        this.regions.add(region);
    }
    private final static int END_5p = 0;
    private final static int END_3p = 1;
    private int[] getGeneCoordinate(){
        int[] coord = new int[2];
      //for(Exon exon: merged){
      //    System.out.println(exon.getChromosome() + ": " + exon.getStart() + " " + exon.getEnd());
      //}
        coord[END_5p] = this.regions.get(0).getStart();
        coord[END_3p] = this.regions.get(this.regions.size() - 1).getEnd();
        return coord;
    }
    private void initRegionConv(){
        int totalLen = 0;
        for(Region region: this.regions){
            totalLen += region.getLength();
        }
        this.totalRegionLength = totalLen;
    }
    // positions only in mapped exons
    private int getMappedRegionX(int logicalPos){
        if(regionMap.containsKey(logicalPos)){
            return regionMap.get(logicalPos);
        }
        throw new RuntimeException("Out of Exon at " + logicalPos);
    }
    private int getRegionX(int logicalPos){
        int offset = 0;
        for(Region region: this.regions){
            if(logicalPos > region.getEnd()){
                offset += region.getLength();
            }else {
                offset += (logicalPos - region.getStart());
                break;
            }
        }
        return (int)(((double)offset/totalRegionLength) * getViewportWidth()) + this.insets.left;
    }
    public void initConv(){
        int[] geneCoord = getGeneCoordinate();
        // int geneLength = endof_3prime - endof_5prime;
        int geneLength = geneCoord[1] - geneCoord[0];
        double[][] matC = {
            {insets.left, insets.top, 1.0},
            {insets.left, this.size.height - insets.bottom, 1.0},
            {this.size.width - insets.right, insets.top, 1.0}
        };
        double[][] matA = {
            {geneCoord[END_5p] - geneLength/9.0, logicalTop, 1},
            {geneCoord[END_5p] - geneLength/9.0, logicalBottom, 1},
            {geneCoord[END_3p]  + geneLength/9.0, 4, 1}
        };
        RealMatrix mA = MatrixUtils.createRealMatrix(matA);
        RealMatrix mC = MatrixUtils.createRealMatrix(matC);
        RealMatrix invA = new LUDecomposition(mA).getSolver().getInverse();
        RealMatrix mB = invA.multiply(mC);

        // dumpMatrix(mB);
        this.convMatrix = mB;

        this.geneCount = this.gene.size();
        this.yoffset = geneCount * (geneThickness + geneMargin);
    }
    private double[] conv(double x, double y){
        double[][] p = {
            {x, y, 1},
            {0, 0, 1},
            {0, 0, 1}
        };
        RealMatrix point = MatrixUtils.createRealMatrix(p);
        RealMatrix result = point.multiply(convMatrix);
        return result.getRow(0);
    }
    private static void dumpMatrix(RealMatrix m) {
        System.out.println("----------------");
        for (int i = 0; i < m.getRowDimension(); i++) {
            System.out.print("{");
            for (int j = 0; j < m.getColumnDimension(); j++) {
                System.out.print(m.getEntry(i, j) + ", ");
            }
            System.out.println("}");
        }
    }
    public void paint(){
        initConv();
        initRegionConv();
        initMapping();
        int height = this.size.height + this.yoffset;
        this.image = new BufferedImage(this.size.width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D)this.image.getGraphics();
        g.setBackground(Color.white);
        g.clearRect(0, 0, this.size.width, height);

        // title
        GenePlotTitle title = new GenePlotTitle(this.sample + " "+ this.gene.get(0).getName());
        System.err.println(this.gene.get(0).getName());
        title.draw(g, this.size.width, height);

        // chromosome
        String chromosome = this.regions.get(0).getChromosome();
        if(!chromosome.toLowerCase().startsWith("chr")){
            chromosome = "chr" + chromosome;
        }
        g.drawString(chromosome, this.insets.left, this.insets.top);

        for(int i = 0; i<this.geneCount; i++){
            RefGene isoform = this.gene.get(i);
        //  for(Exon exon: isoform.getExons()){
        //      System.out.println("exon: "+ exon.getStart() + " "  + exon.getEnd());
        //  }
            drawTrackLabel(g, isoform, i);
            drawIntron(g, isoform, i);
            drawExon(g, isoform, i);
        }
        drawScoreAxes(g);
        drawGenomeAxis(g);
        drawMapping(g);
        drawRegions(g);
        drawCNVs(g);

        // drawScoreTest(g);
    }
    private void drawScoreTest(Graphics2D g){
        int pos = 89350000;
        double score = 0.7;
        drawScore(g, pos, score);
        
    }
    private void drawRegions(Graphics2D g){
        for(Region region: this.regions){
            double[] ratio = region.getRatioNormalized(sample);
            /*
            if(region.getStart() == 157222343){
                for(int i = 0; i<ratio.length; i++){
                    System.out.println("doubt: " + ratio[i]);
                }
            }*/
            for(int pos = region.getStart(); pos < region.getEnd(); pos++){
                try {
                    drawScore(g, pos, ratio[pos - region.getStart()]);   
                }catch(Exception e){
                    // e.printStackTrace(); ignore NullPointerException
                }
            }
        }
    }
    private void drawCNVs(Graphics2D g){
        for(CNV cnv: this.cnvs){
            double[] ratio = cnv.getRatioNormalized(sample, cnv.getStart(), cnv.getEnd());
            for(int pos = cnv.getStart(); pos < cnv.getEnd(); pos++){
                try {
                    drawCNVScore(g, pos, ratio[pos - cnv.getStart()], colorMap.get(cnv.getType()));   
                }catch(Exception e){
                    // e.printStackTrace(); ignore NullPointerException
                }
            }
        }
    }
    private void drawCNVScore(Graphics2D g, int pos, double score, Color color){
        double[] p = conv(0, score);
        p[0] = getMappedRegionX(pos);
        final double size = 4.0;
        // System.out.println("point = " + p[0] + " "+ p[1]);
        java.awt.geom.Ellipse2D.Double shape = new java.awt.geom.Ellipse2D.Double(p[0] - size/2, p[1] - size/2 + this.yoffset, size, size);
        g.setColor(color);
        g.draw(shape);
        g.setColor(Color.black);
    }
    private void drawScore(Graphics2D g, int pos, double score){
        double[] p = conv(0, score);
        p[0] = getMappedRegionX(pos);
        final double size = 4.0;
        // System.out.println("point = " + p[0] + " "+ p[1]);
        java.awt.geom.Ellipse2D.Double shape = new java.awt.geom.Ellipse2D.Double(p[0] - size/2, p[1] - size/2 + this.yoffset, size, size);
        g.draw(shape);
    }
    private void initMapping(){
        int offset = 0;
        for(Region region: this.regions){
            for(int pos = region.getStart(); pos < region.getEnd(); pos++){
                int phys_pos = (int)(((double)(offset+pos-region.getStart())/totalRegionLength) * getViewportWidth()) + this.insets.left;
                regionMap.put(pos, phys_pos);
            }
            offset += region.getLength();
        }
    }
    private int getViewportWidth(){
        return this.size.width - this.insets.right - this.insets.left;
    }
    private void drawMapping(Graphics2D g){
        for(Region region: this.regions){
            // System.err.println(region.getStart() + " " + region.getEnd());
            double[] p1 = conv(region.getStart(), REFGENE_TRACK_BASE_POS);
            p1[1] += this.yoffset;

            final double ypos = 2.1;
            double p2X = getRegionX(region.getStart());
            double p2Y = conv(0, ypos)[1];
            g.setColor(Color.gray);
            Shape line = new Line2D.Double(p1[0], p1[1], p2X, p2Y + this.yoffset);
            g.draw(line);
        }
    }
    private void drawGenomeAxis(Graphics2D g){
        int[] geneCoord = getGeneCoordinate();
        int geneLength = geneCoord[1] - geneCoord[0];
        // main axis
        double[] p1 = conv(geneCoord[END_5p] - (geneLength/9.0), GENOME_AXIS_POS);
        double[] p2 = conv(geneCoord[END_3p] + (geneLength/9.0), 0);
        GenePlotGenomeAxis axis = new GenePlotGenomeAxis();
        // axis.draw(g, (int)p1[0], (int)p2[0], this.insets.left, (int)p1[1], (int)(this.size.width - this.insets.left - this.insets.right));
        axis.draw(g, geneCoord[END_5p] - (int)(geneLength/9.0), geneCoord[END_3p] + (int)(geneLength/9.0), this.insets.left, (int)p1[1], (int)(this.size.width - this.insets.left - this.insets.right));
    }
    private void drawScoreAxes(Graphics2D g){
        // use this.yoffset
        int[] geneCoord = getGeneCoordinate();
        // main axis
        double[] p1_y0 = conv(geneCoord[END_5p], 0);
        double[] p2_y0 = conv(geneCoord[END_3p], 0);
        double[] p1_yp1 = conv(geneCoord[END_5p], 1);
        double[] p2_yp1 = conv(geneCoord[END_3p], 1);
        double[] p1_ym1 = conv(geneCoord[END_5p], -1);
        double[] p2_ym1 = conv(geneCoord[END_3p], -1);
        Shape line0 = new Line2D.Double(this.insets.left, p1_y0[1] + this.yoffset, this.size.width - this.insets.right, p2_y0[1] + this.yoffset);
        Shape line_plus1 = new Line2D.Double(this.insets.left, p1_yp1[1] + this.yoffset, this.size.width - this.insets.right, p2_yp1[1] + this.yoffset);
        Shape line_minus1 = new Line2D.Double(this.insets.left, p1_ym1[1] + this.yoffset, this.size.width - this.insets.right, p2_ym1[1] + this.yoffset);
        g.setColor(Color.black);
        g.draw(line0);
        Stroke defaultStroke = g.getStroke();
        Stroke dashedStroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{9}, 0);
        g.setColor(Color.gray);
        g.setStroke(dashedStroke);
        g.draw(line_plus1);
        g.draw(line_minus1);
        g.setStroke(defaultStroke);
        g.setColor(Color.black);

        double[] p1_yp2 = conv(geneCoord[END_5p], 2);
        double[] p1_ym2 = conv(geneCoord[END_5p], -2);
        Shape yaxis = new Line2D.Double(this.insets.left, p1_yp2[1] + this.yoffset, this.insets.left, p1_ym2[1] + this.yoffset);
        g.draw(yaxis);

        java.awt.FontMetrics fm = g.getFontMetrics();
        String tick0 = "0.0";
        String tickp1 = "1.0";
        String tickm1 = "-1.0";
        String tickp2 = "2.0";
        String tickm2 = "-2.0";

        int stringHeight = fm.getHeight();
        g.drawString(tick0, (float)insets.left - fm.stringWidth(tick0) - 2, (float)p2_y0[1] + this.yoffset + stringHeight/2);
        g.drawString(tickp1, (float)insets.left - fm.stringWidth(tickp1) - 2, (float)p1_yp1[1] + this.yoffset + stringHeight/2);
        g.drawString(tickp2, (float)insets.left - fm.stringWidth(tickp2) - 2, (float)p1_yp2[1] + this.yoffset + stringHeight/2);
        g.drawString(tickm1, (float)insets.left - fm.stringWidth(tickm1) - 2, (float)p1_ym1[1] + this.yoffset + stringHeight/2);
        g.drawString(tickm2, (float)insets.left - fm.stringWidth(tickm2) - 2, (float)p1_ym2[1] + this.yoffset + stringHeight/2);
        // tick at 2.0
        Shape tick = new Line2D.Double(this.insets.left, p1_yp2[1] + this.yoffset, this.insets.left + 2, p1_yp2[1] + this.yoffset);
        g.draw(tick);

        // tick at -2.0
        tick = new Line2D.Double(this.insets.left, p1_ym2[1] + this.yoffset, this.insets.left + 2, p1_ym2[1] + this.yoffset);
        g.draw(tick);
    }
    private void drawTrackLabel(Graphics2D g, RefGene isoform, int index){
        double offset = index * (geneMargin + geneThickness);
        java.awt.FontMetrics fm = g.getFontMetrics();
        double x = this.insets.left  - fm.stringWidth(isoform.accession) - 4;
        double y = offset + fm.getHeight()/2.0;
        double[] p = conv(0, REFGENE_TRACK_BASE_POS); // 0 is dummy
        g.setColor(Color.black);
        g.drawString(isoform.accession,(float) x, (float)(p[1] + y) - fm.getDescent());
    }
    private void drawExon(Graphics2D g, RefGene isoform, int index){
        double offset = index * (geneMargin + geneThickness);
        for(Exon exon: isoform.getExons()){
            double[] p1 = conv(exon.getStart(), REFGENE_TRACK_BASE_POS);
            double[] p2 = conv(exon.getEnd(), REFGENE_TRACK_BASE_POS);
            g.setColor(geneBackgroundColor);
            Shape exonShape = new Rectangle2D.Double(p1[0], p1[1] - geneThickness/2.0 + offset, p2[0] - p1[0], geneThickness);
            g.fill(exonShape);
            exonShape = new Rectangle2D.Double(p1[0], p1[1] - geneThickness/2.0 + offset, p2[0] - p1[0], geneThickness);
            g.setColor(geneBorderColor);
            g.draw(exonShape);

        }
    }
    private void drawIntron(Graphics2D g, RefGene isoform, int index){
        double offset = index * (geneMargin + geneThickness);
        double[] p1 = conv(isoform.getFirstExon().getStart(), REFGENE_TRACK_BASE_POS);
        double[] p2 = conv(isoform.getLastExon().getEnd(), REFGENE_TRACK_BASE_POS);
        // System.out.println("src geom: " + isoform.getFirstExon().getStart()  + " "+ isoform.getLastExon().getEnd());
        // System.out.println("dst geom: " + p1[0] + " " + p1[1] + " " + p2[0] + " " + p2[1]);

        Shape line = new Line2D.Double(p1[0], p1[1] + offset, p2[0], p2[1] + offset);
        // System.out.println("geom: " + p1[0] + " " + (p1[1] + offset) + " " + p2[0] + " " + (p2[1] + offset));
        g.setColor(intronColor);
        g.draw(line);
    }
    public boolean write(String filename, String format) throws java.io.IOException {
        return javax.imageio.ImageIO.write(this.image, format, new java.io.FileOutputStream(filename));
    }

    public static void main(String[] argv){
        try {
            String testGene = "ANKRD11";
            ArrayList<RefGene> genes = RefGene.load(new java.util.zip.GZIPInputStream(new java.io.FileInputStream("./resources/UCSC-20181203/refGene.txt.gz")));
            java.util.HashMap<String, ArrayList<RefGene>> map = new java.util.HashMap<>();
            for(RefGene gene: genes){
                if(!map.containsKey(gene.getName())){
                    map.put(gene.getName(), new ArrayList<RefGene>());
                }
                map.get(gene.getName()).add(gene);
            }
            ArrayList<RefGene> targetGene = map.get(testGene);
            // System.out.println(testGene + " has " + targetGene.size() + " isoforms");
            
            // for the sake of test
            GenePlot plot = new GenePlot("DUMMY SAMPLE", targetGene, 800, 600);
            plot.paint();
            plot.write("test.png", "png");
            /*
            ArrayList<String> genes = new ArrayList<String>();
            genes.add("ANKRD11");
            File file = new File("resources/S04380110_Padded.bed");
            SureSelect ss = new SureSelect(file, genes, true);
            for(String geneName: genes){
                Gene gene = new Gene(geneName, ss.getRegions(geneName));
            }
            */
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
