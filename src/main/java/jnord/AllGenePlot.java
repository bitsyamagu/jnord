package jnord;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.awt.Color;
import java.awt.geom.Dimension2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.data.DataSeries;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.points.PointRenderer;
import de.erichseifert.gral.plots.points.DefaultPointRenderer2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.plots.axes.Axis;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.graphics.Insets2D;
// import de.erichseifert.gral.graphics.Label;
import de.erichseifert.gral.graphics.Orientation;
import java.awt.geom.Rectangle2D;
import java.awt.Stroke;
import java.awt.BasicStroke;


public class AllGenePlot {
    String sample = null;
    boolean showLegend = true;
    static Color[] colors = null;

    static {
        System.setProperty("java.awt.headless", "true");
    };
    public AllGenePlot(String sampleName, boolean showLegend){
        this.sample = sampleName;
        this.showLegend = showLegend;
    }
    static Rectangle2D.Double rect = new Rectangle2D.Double(0.0, 0.0, 4.0, 4.0);
    static Stroke stroke = new BasicStroke(4.0f);
    public void plot(ArrayList<Gene> genes) throws FileNotFoundException, IOException{
        int pos = 0;
        // ArrayList<DataTable> allData = new ArrayList<>();
        XYPlot plot = new XYPlot();
        plot.setBackground(Color.WHITE);
        plot.setInsets(new Insets2D.Double(20.0, 40.0, 40.0, 20.0));
        plot.getTitle().setText(sample + " All Genes");
        if(showLegend){
            plot.setLegendVisible(true);
            plot.getLegend().setOrientation(Orientation.HORIZONTAL);
        }
        
        IwamaColorScheme colorScheme = new IwamaColorScheme();
        if(colors == null){
            colors = new Color[genes.size()];
            if(AnalysisContext.randomColorScheme){
                ColorScheme.setSaturation(200);
            }
            for(int i = 0; i<colors.length; i++){
                if(AnalysisContext.randomColorScheme){
                    colors[i] = ColorScheme.randomColor();
                }else {
                    colors[i] = colorScheme.nextColor();
                }
            }
        }
        int colorIndex = 0;
        
        // PrintWriter pw = new PrintWriter(new FileWriter(this.sample +".summary.txt"));
        for(Gene gene: genes){
            // Label label = new Label(gene.getName());
            // Dimension2D labelSize = label.getPreferredSize();
            // label.setBounds(pos, 3.0, labelSize.getWidth(), labelSize.getHeight());
            DataTable geneDataTable = new DataTable(Integer.class, Double.class);
            for(Region region: gene.getRegions()){
                // pw.println(region.getChromosome() + "\t" + region.getStart() + "\t" + region.getEnd() + "\t" + gene.getName());
                double[] ratio = region.getRatioNormalized(this.sample);
                double[] sn = region.getNormalizedSN();
             // for(double r: ratio){
             //     geneDataTable.add(pos++, r);
             // }
                for(int i = 0; i<ratio.length; i++){
                    if(sn[i] > AnalysisContext.SN_threshold){
                        geneDataTable.add(pos++, ratio[i]);
                    }
                }
            }
            DataSeries series = null;
            plot.add(series = new DataSeries(gene.getName(), geneDataTable));
            // plot.add(label);
            LineRenderer line_renderer = new DefaultLineRenderer2D();
            PointRenderer point_renderer = new DefaultPointRenderer2D();
            point_renderer.setShape(rect);
            line_renderer.setStroke(stroke);
            line_renderer.setColor(colors[colorIndex]);
            point_renderer.setColor(colors[colorIndex]);
            point_renderer.setShape(new java.awt.geom.Line2D.Double());
            plot.setLineRenderers(series, line_renderer);
            plot.setPointRenderers(series, point_renderer);
            // allData.add(geneDataTable);
            colorIndex++;
        }
        Axis a = plot.getAxis(XYPlot.AXIS_Y);
        a.setMax(4.0);
        a.setMin(0.0);
        // pw.close();
        DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
        writer.write(plot, new FileOutputStream(new File(sample  + ".all_genes.png")), 1900, 600);

        PrintWriter pw  = new PrintWriter(new BufferedWriter(new FileWriter(sample + ".all_genes.txt")));
        for(Gene gene: genes){
            Region firstExon = gene.getFirst();
            Region lastExon = gene.getLast();
            pw.println(gene.getName()+ "\t"  + firstExon.getChromosome() + "\t" + firstExon.getStart() + "\t" + lastExon.getEnd());
        }
        pw.close();
        System.err.println("gene info are written in " + sample + ".all_genes.txt");
    }
}
