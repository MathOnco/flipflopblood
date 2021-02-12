package Epigenome;

//to make another EpigenomeSkeleton
public class IntArrayEpigenome extends EpigenomeSkeleton {
    short[][]data;

    public IntArrayEpigenome(int chromNumber, int length) {
        data=new short[chromNumber][length];
    }

    @Override
    public boolean Get(int chrom, int i) {
        return data[chrom][i] != 0;
    }

    @Override
    public void Set(int chrom, int i, boolean val) {
        data[chrom][i]=val?(short)1:(short)0;
    }

    @Override
    public EpigenomeSkeleton Copy() {
        IntArrayEpigenome ret=new IntArrayEpigenome(data.length, data[0].length);
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i],0,ret.data[i],0,data[i].length);
        }
        return ret;
    }

    @Override
    public int numChroms() {
        return data.length;
    }

    @Override
    public int lengthChrom(int i) {
        return data[0].length;
    }

}
