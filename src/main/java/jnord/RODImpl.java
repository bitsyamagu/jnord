package jnord;

import java.util.HashMap;

public class RODImpl implements ROD, Comparable<ROD>{
    String chr;
    int start;
    int end;
    public static final HashMap<String, Integer> ROD_INDEX = new HashMap<String, Integer>();
    static {
        ROD_INDEX.put("1", 1);
        ROD_INDEX.put("2", 2);
        ROD_INDEX.put("3", 3);
        ROD_INDEX.put("4", 4);
        ROD_INDEX.put("5", 5);
        ROD_INDEX.put("6", 6);
        ROD_INDEX.put("7", 7);
        ROD_INDEX.put("8", 8);
        ROD_INDEX.put("9", 9);
        ROD_INDEX.put("10", 10);
        ROD_INDEX.put("11", 11);
        ROD_INDEX.put("12", 12);
        ROD_INDEX.put("13", 13);
        ROD_INDEX.put("14", 14);
        ROD_INDEX.put("15", 15);
        ROD_INDEX.put("16", 16);
        ROD_INDEX.put("17", 17);
        ROD_INDEX.put("18", 18);
        ROD_INDEX.put("19", 19);
        ROD_INDEX.put("20", 20);
        ROD_INDEX.put("21", 21);
        ROD_INDEX.put("22", 22);
        ROD_INDEX.put("X", 23);
        ROD_INDEX.put("Y", 24);
    };
    public RODImpl(String chr, int start, int end){
        this.chr = chr;
        this.start = start;
        this.end = end;
    }
    public void setChromosome(String c){
        this.chr = c;
    }
    public void setStart(int s){
        this.start = s;
    }
    public void setEnd(int e){
        this.end = e;
    }
    public String getChromosome(){
        return this.chr;
    }
    public int getStart(){
        return this.start;
    }
    public int getEnd(){
        return this.end;
    }
    public int getLength(){
        return this.end - this.start + 1;
    }
    public boolean hasIntersection(ROD rod){
        ROD upstream = null;
        ROD downstream = null;
        if(this.chr.equals(rod.getChromosome())){
            if(this.start < rod.getStart() || (this.start == rod.getStart() && this.end < rod.getEnd())){
                upstream = this;
                downstream = rod;
            }else {
                upstream = rod;
                downstream = this;
            }
            if(upstream.getEnd() > downstream.getStart() && upstream.getEnd() <= downstream.getEnd()){
                return true;
            }else if(downstream.getEnd() > upstream.getStart() && downstream.getEnd() <= upstream.getEnd()){
                return true;
            }
        }
        return false;
    }
    public int compareTo(ROD rod){
        Integer aindex = ROD_INDEX.get(chr);
        Integer bindex = ROD_INDEX.get(rod.getChromosome());
        if(aindex != null && bindex != null){
            int c = aindex.compareTo(bindex);
            if(c != 0){
                return c;
            }
            if(start < rod.getStart()){
                return -1;
            }else if(start > rod.getStart()){
                return 1;
            }
            if(end < rod.getEnd()){
                return -1;
            }else if(end > rod.getEnd()){
                return 1;
            }
        }else if(aindex != null && bindex == null){
            return -1;
        }else if(aindex == null && bindex != null){
            return 1;
        }
        return 0;
    }
    public void merge(ROD rod){
        if(this.start > rod.getStart()){
            this.start = rod.getStart();
        }
        if(this.end < rod.getEnd()){
            this.end = rod.getEnd();
        }
    }
    public boolean equals(ROD rod){
        if(chr.equals(rod.getChromosome()) && start == rod.getStart() && end == rod.getEnd()){
            return true;
        }
        return false;
    }
}
