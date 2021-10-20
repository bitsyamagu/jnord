package jnord;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.ReferenceSequence;
import java.io.File;
import java.io.FileNotFoundException;

public class GCPercent {
    IndexedFastaSequenceFile fasta;
    public GCPercent(File file, File indexFile) throws FileNotFoundException{
        fasta = new IndexedFastaSequenceFile(file, new FastaSequenceIndex(indexFile));
    }
    public double[] getGCpercent100(Region r){
        double[] buf = new double[r.getEnd() - r.getStart() + 1];
        ReferenceSequence ref = fasta.getSubsequenceAt(r.getChromosome(), r.getStart()-50, r.getEnd()+50+1);
        byte[] bases = ref.getBases();
        // System.err.println("loading gc percent for " + r.getChromosome() + " " + r.getStart() + " " + r.getEnd());
        for(int i = 0; i< r.getEnd() - r.getStart() + 1; i++){
            int count = 0;
            for(int j = 0; j<100; j++){
                if(bases[i+j] == 'G' || bases[i+j] == 'C' || bases[i+j] == 'g' || bases[i+j] == 'c'){
                    count++;
                }
            }
            buf[i] = count/100.0f;
        }
        // dump(r, buf);
        return buf;
    }
    public void dump(Region r, double[] buf){
        System.err.println("loaded gc percent for " + r.getChromosome() + " " + r.getStart() + " " + r.getEnd());
        for(int i = 0; i<buf.length; i++){
            System.err.println(r.getChromosome() + " " + (r.getStart() + i) + " " + buf[i]);
        }
    }
}
