package jnord;

public interface ROD {
    public String getChromosome();
    public int getStart();
    public int getEnd();
    public void merge(ROD rod);
}
