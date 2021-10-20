package jnord;

import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Color;

public class GenePlotTitle {
    String titleText = "";
    public GenePlotTitle(String text){
        this.titleText = text;
    }
    public String getText(){
        return this.titleText;
    }
    public void draw(Graphics2D g2d, int width, int height){
        int centerX = width/2;
        FontMetrics fm = g2d.getFontMetrics();
        int strHeight = fm.getHeight();
        int strWidth = fm.stringWidth(this.titleText);

        g2d.setColor(Color.black);
        g2d.drawString(this.titleText, centerX - strWidth/2, strHeight + 2);
    }
}
