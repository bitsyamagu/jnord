package jnord;

import java.awt.Color;
import java.util.ArrayList;

public class IwamaColorScheme {
    static int saturation = 255;
    static ArrayList<Color> colors = new ArrayList<>();
    int index = -1;
    static {
        colors.add(new Color(0xff, 0x00, 0x00));
        colors.add(new Color(0x7f, 0xbf, 0xff));
        colors.add(new Color(0xff, 0xff, 0x00));
        colors.add(new Color(0xbf, 0x7f, 0xff));
        colors.add(new Color(0x00, 0xff, 0x00));
        colors.add(new Color(0xff, 0x7f, 0xbf));
        colors.add(new Color(0x00, 0xff, 0xff));
        colors.add(new Color(0xff, 0xbf, 0x7f));
        colors.add(new Color(0x00, 0x00, 0xff));
        colors.add(new Color(0xbf, 0xff, 0x7f));

        colors.add(new Color(0xff, 0x00, 0xff));
        colors.add(new Color(0x7f, 0xff, 0xbf));
        colors.add(new Color(0xff, 0x7f, 0x00));
        colors.add(new Color(0x7f, 0x7f, 0xff));
        colors.add(new Color(0x7f, 0xff, 0x00));
        colors.add(new Color(0xff, 0x7f, 0xff));
        colors.add(new Color(0x00, 0xff, 0x7f));
        colors.add(new Color(0xff, 0x7f, 0x7f));
        colors.add(new Color(0x00, 0x7f, 0xff));
        colors.add(new Color(0xff, 0xff, 0x7f));

        colors.add(new Color(0x7f, 0x00, 0xff));
        colors.add(new Color(0x7f, 0xff, 0x7f));
        colors.add(new Color(0xff, 0x00, 0x7f));
        colors.add(new Color(0x7f, 0xff, 0xff));
        colors.add(new Color(0x60, 0x60, 0x60));
    };
    public Color nextColor(){
        index++;
        return colors.get(index%colors.size());
    }
}
