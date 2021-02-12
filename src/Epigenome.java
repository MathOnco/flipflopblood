import Epigenome.EpigenomeSkeleton;
import Epigenome.GenomeInfo;
import HAL.Rand;

public class Epigenome extends GenomeInfo<Epigenome> {
    float h, s, v;
    EpigenomeSkeleton epiData;
    static double[][] EPICHANGERATE;
    Rand rn;
    boolean mutated;
    static int[][] binsites;

    Epigenome(float h, float s, float v, EpigenomeSkeleton ParentEpiGenomeArray, double[][] EPICHANGERATE, Rand rn, int[][] BinnedEpiSites) {
        this.h = h;
        this.s = s;
        this.v = v;
        this.epiData = ParentEpiGenomeArray;
        this.EPICHANGERATE = EPICHANGERATE;
        this.rn = rn;
        this.binsites = BinnedEpiSites;
    }

    @Override
    public Epigenome _RunPossibleMutation(){
        int triggerTrack=0;
        EpigenomeSkeleton newArray= epiData.Copy();

        for (int i = 0; i < epiData.numChroms(); i++) {
            for (int j = 0; j < epiData.lengthChrom(i); j++) {
                if(newArray.Get(i,j)==true){ // CpG site is on
                    if(rn.Double()<=this.EPICHANGERATE[1][binsites[0][j]]){
                        triggerTrack++;
                        newArray.Flip(i,j);
                    }
                } else {
                    if(rn.Double()<=this.EPICHANGERATE[0][binsites[0][j]]){
                        triggerTrack++;
                        newArray.Flip(i,j);
                    }
                }
            }
        }
        if(triggerTrack>0){
            mutated = true;
            return new Epigenome((float)rn.Double(), 1f , 0.75f, newArray, this.EPICHANGERATE, rn, binsites);
        } else {
            mutated = false;
            return null;
        }
    }

    /*
    Don't use
     */
    public String GenomeInfoStr() {
        return null;
    }
}