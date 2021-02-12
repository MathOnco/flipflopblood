package Epigenome;

import HAL.Rand;

public abstract class EpigenomeSkeleton {

    //Get and Set are "atomic" functions required for the interface to work

    public abstract boolean Get(int chrom, int i);

    public abstract void Set(int chrom, int i, boolean val);

    public abstract EpigenomeSkeleton Copy();

    public abstract int numChroms();
    public abstract int lengthChrom(int i);

    //default functions use "atomic" functions to do things, and don't have to be implemented

    public int GetInt(int chrom, int i){
        return Get(chrom, i) ? 1:0;
    }

    public void Flip(int chrom, int i){
        Set(chrom, i,!Get(chrom, i));
    }

    public EpigenomeSkeleton InitializeMethylationState(Rand rn){
        for (int chrom = 0; chrom < numChroms(); chrom++) {
            for (int site = 0; site < lengthChrom(chrom); site++) {
                Set(chrom,site,false);
            }
        }
        return this;
    }

}
