package jnord;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.awt.Color;
import java.util.ArrayList;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.points.PointRenderer;
import de.erichseifert.gral.plots.points.DefaultPointRenderer2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.graphics.Insets2D;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class InvariantSetNormalizationPlot {
    String sample = null;

    static {
        System.setProperty("java.awt.headless", "true");
    };

    public InvariantSetNormalizationPlot(String sampleName){
        this.sample = sampleName;
    }
    public void plot(ArrayList<Gene> genes, PolynomialSplineFunction curve) throws IOException, FileNotFoundException{
        File preNorm = new File(this.sample + "_before_invnorm.png");
        File postNorm = new File(this.sample + "_after_invnorm.png");
        DataTable preTable = new DataTable(Double.class, Double.class);
        DataTable postTable = new DataTable(Double.class, Double.class);

        for(Gene gene: genes){
            for(Region region: gene.getRegions()){
                double[] medians = region.getMedian();
                double[] pre = region.getDepth(sample);
                double[] post = region.getNormalizedBySample(sample);
                for(int i = 0; i<medians.length; i++){
                    preTable.add(Math.log(medians[i]), Math.log(pre[i]));
                    postTable.add(Math.log(medians[i]), Math.log(post[i]));
                }
            }
        }
        XYPlot plot = new XYPlot(preTable);
        plot.getTitle().setText(sample + " before normalization");
        plot.setInsets(new Insets2D.Double(20.0, 40.0, 40.0, 20.0));
        plot.setBackground(Color.WHITE);

        double[] knots = curve.getKnots();
        DataTable splineData = new DataTable(Double.class, Double.class);
        for(int i = 0; i<knots.length; i++){
            splineData.add(knots[i], curve.value(knots[i]));
        }
        LineRenderer lineRenderer = new DefaultLineRenderer2D();
        PointRenderer pointRenderer = new DefaultPointRenderer2D();
        lineRenderer.setColor(Color.red);
        pointRenderer.setColor(Color.red);
        plot.add(splineData);
        plot.setLineRenderers(splineData, lineRenderer);
        plot.setPointRenderers(splineData, pointRenderer);

        DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
        // writer.write(plot, new FileOutputStream(preNorm), 800, 600);
        writer.write(plot, new FileOutputStream(preNorm), 1900, 600);

        plot = new XYPlot(postTable);
        plot.getTitle().setText(sample + " after normalization");
        plot.setInsets(new Insets2D.Double(20.0, 40.0, 40.0, 20.0));
        plot.setBackground(Color.WHITE);
        writer = DrawableWriterFactory.getInstance().get("image/png");
        // writer.write(plot, new FileOutputStream(postNorm), 800, 600);
        writer.write(plot, new FileOutputStream(postNorm), 1900, 600);
    }
    
}
