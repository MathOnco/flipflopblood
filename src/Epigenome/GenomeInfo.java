package Epigenome;

/**
 * should be declared myType extends GenomeInfo <myType>
 */
public abstract class GenomeInfo <T extends GenomeInfo> {
    public int id;
    int popSize;
    T next;
    T prev;
    GenomeTracker myTracker;


    /**
     * gets the current number of clones that share this genome
     */
    public int GetClonePop(){
        return popSize;
    }

    /**
     * ignore
     */
    void _Init(GenomeTracker myTracker, int id, T next, T prev){
        this.myTracker=myTracker;
        this.id=id;
        this.next=next;
        this.prev=prev;
        this.popSize=0;
    }

    /**
     * removes clone from GenomeInfo population
     */
    public void DisposeClone(){
        myTracker.DisposeClone(this);
    }

    /**
     * adds new clone to GenomeInfo population
     */
    public T NewChild(){
        popSize++;
        return (T)this;
    }

    /**
     * returns a GenomeInfo instance that is either identical to the original or a mutant
     */
    public T PossiblyMutate(){
        T nextGenome= _RunPossibleMutation();
        if(nextGenome==null){
            return (T)this;
        }
        DisposeClone();
        myTracker.AddMutant(this,nextGenome);
        return nextGenome;
    }
    public String FullLineageInfoStr(String delim){
        return myTracker.FullLineageInfoStr(id,delim);
    }

    /**
     * a potential mutation event, return null if the genome did not change, otherwise return a new genome with the change inside
     * do not change the calling genome!!!!!!!!!!!!!!
     */
    public abstract T _RunPossibleMutation();

    /**
     * returns a string with info about the genome to be stored
     */
    public abstract String GenomeInfoStr();

    public int IDGetter(){ return id; }
}
