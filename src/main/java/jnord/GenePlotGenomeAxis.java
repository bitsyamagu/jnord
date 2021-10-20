package jnord;

// import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Color;
import java.awt.geom.Line2D;
import java.awt.Shape;


public class GenePlotGenomeAxis {
    public GenePlotGenomeAxis(){
    }
    public void draw(Graphics2D g2d, int logicalStart, int logicalEnd, int x, int y, int length){
        // FontMetrics fm = g2d.getFontMetrics();
        // int strHeight = fm.getHeight();
        // int strWidth = fm.stringWidth(this.titleText);

        g2d.setColor(Color.black);
        Shape line = new Line2D.Double(x, y, x + length, y);
        g2d.draw(line);

        // 223k
        int logicalLen = logicalEnd - logicalStart;
        // g2d.drawString(this.titleText, centerX - strWidth/2, strHeight + 2);
        // 22k
        final int[] ranges = {10, 100, 1000, 10000, 100000, 1000000, 10000000};
        int majorTickInterval = ranges[(int)Math.log10(logicalLen)-1];
        int minorTickInterval = majorTickInterval/10;
        // System.err.println("majorTickInterval: " + majorTickInterval);
        // System.err.println("logical length: " + (logicalEnd - logicalStart));
        int majorTickStart = 0;
        int minorTickStart = 0;
        for(int i = logicalStart; i<logicalEnd; i++){
            if((i % majorTickInterval) == 0){
                majorTickStart = i;
                break;
            }
        }
        for(int i = logicalStart; i<logicalEnd; i++){
            if((i % minorTickInterval) == 0){
                minorTickStart = i;
                break;
            }
        }
        // System.out.println("majorTickStart " + majorTickStart);
        java.awt.FontMetrics fm = g2d.getFontMetrics();
        for(int majorTick = majorTickStart; majorTick < logicalEnd; majorTick += majorTickInterval){
            double xmtick = ((double)(majorTick - logicalStart))/(logicalEnd - logicalStart) * (double)length + x;
            // label
            // String label = String.valueOf(majorTick);
            String label = java.text.NumberFormat.getInstance().format(majorTick);
            int labelWidth = fm.stringWidth(label);
            int labelPos = (int)xmtick - labelWidth/2;
            g2d.drawString(label, labelPos, y - fm.getHeight() - 2);

            // System.out.println(xmtick + " " + y);
            line = new Line2D.Double(xmtick, y, xmtick, y+8);
            g2d.draw(line);
        }
        for(int minorTick = minorTickStart; minorTick < logicalEnd; minorTick += minorTickInterval){
            double xmtick = ((double)(minorTick - logicalStart))/(logicalEnd - logicalStart) * (double)length + x;
            line = new Line2D.Double(xmtick, y, xmtick, y+2);
            g2d.draw(line);
        }
    }
}
