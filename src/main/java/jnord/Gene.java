package jnord;

import java.util.Collection;
import java.util.ArrayList;

public class Gene implements Comparable<Gene> {
    String name = null;
    ArrayList<Region> regions = new ArrayList<>();
    
    public Gene(String name){
        this.name = name;
    }
    public Gene(String name, java.util.List<Region> list){
        this.name = name;
        regions.addAll(list);
     // for(Region r: list){
     //     regions.add(r);
     // }
        java.util.Collections.sort(regions);
        // System.out.println(" " + regions.size() );
    }
    public String getName(){
        return this.name;
    }
    public Region getFirst(){
        return regions.get(0);
    }
    public Region getLast(){
        return regions.get(regions.size() -1);
    }
    public java.util.List<Region> getRegions(){
        return this.regions;
    }
    public int compareTo(Gene another){
        return regions.get(0).getStart() - another.regions.get(0).getStart();
    }
}
