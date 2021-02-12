package Epigenome;

import HAL.Tools.FileIO;

import java.util.ArrayList;

/**
 * should be instantiated GenomeTracker<myGenomeInfoType>
 */
public class GenomeTracker <T extends GenomeInfo>{
    ArrayList<Integer> parentIDs;
    ArrayList<String> allGenomeInfos;
    ArrayList<int[]> CloneCounts;
    T livingCloneInfos;
    final T progenitors;
    int nMutations=0;
    int nActiveClones=0;

    /**
     * creates a new genome tracker, with a progenitor genome that will be returned whenever NewProgenitor is called
     * @param progenitorGenome the starting genome
     * @param trackParents whether to track the lineages
     * @param trackGenomeInfos whether to track the genome information
     */
    public GenomeTracker(T progenitorGenome,boolean trackParents,boolean trackGenomeInfos){
        CloneCounts=new ArrayList<>();
        progenitors = progenitorGenome;
        if(trackParents) {
            parentIDs = new ArrayList<>();
            parentIDs.add(-1);
        }
        if(trackGenomeInfos) {
            allGenomeInfos = new ArrayList<>();
            allGenomeInfos.add(progenitors.GenomeInfoStr());
        }
        progenitorGenome._Init(this,0,null,null);
    }

    /**
     * returns the progenitor genome, and increases the population size
     */
    public T NewProgenitor(){
        if(progenitors.popSize==0){
            nActiveClones++;
            if(livingCloneInfos!=null) { livingCloneInfos.prev = progenitors; }
            progenitors.next=livingCloneInfos;
            progenitors.prev=null;
            livingCloneInfos =progenitors;
        }
        progenitors.popSize++;
        return progenitors;
    }

    /**
     * ignore
     */
    T AddMutant(T parent,T child) {
        nMutations++;
        nActiveClones++;
        if(allGenomeInfos!=null) { allGenomeInfos.add(child.GenomeInfoStr()); }
        if(parentIDs!=null) { parentIDs.add(parent.id); }
        child._Init(this,NumMutations(), livingCloneInfos,null);
        child.popSize++;
        if(livingCloneInfos!=null) { livingCloneInfos.prev = child; }
        livingCloneInfos =child;
        return child;
    }

    /**
     * returns the total number of unique clone populations that have ever existed
     */
    public int NumMutations(){
        return nMutations;
    }

    /**
     * returns the total number of unique clone populations that currently exist
     */
    public int NumActiveMutantTypes(){
        return nActiveClones;
    }

    void DisposeClone(T cloneInfo){
        cloneInfo.popSize--;
        if(cloneInfo.popSize<0){
            System.out.println("Error! a clone popsize is less than 0!");
        }
        if(cloneInfo.popSize==0){
            nActiveClones--;
            //remove livingCloneInfos with 0 population
            if(livingCloneInfos ==cloneInfo){
                if(livingCloneInfos.prev!=null){
                    System.out.println("something here!");
                }
                livingCloneInfos =(T)cloneInfo.next;
            }
            if(cloneInfo.next!=null){
                cloneInfo.next.prev=cloneInfo.prev;
            }
            if(cloneInfo.prev!=null){
                cloneInfo.prev.next=cloneInfo.next;
            }
        }
    }

    /**
     * records the sizes of all clonal populations that are alive when the function is called
     */
    public void RecordClonePops(){
        int[]ret=new int[nActiveClones*2];
        int iClone=0;
//        if(progenitors.popSize!=0){
//            ret[0]=progenitors.id;
//            ret[1]=progenitors.popSize;
//            iClone++;
//        }
        T clone= livingCloneInfos;
        if(livingCloneInfos!=null && livingCloneInfos.prev!=null){
            System.out.println("here?");
        }
        while(clone !=null){
            ret[iClone*2]=clone.id;
            ret[iClone*2+1]=clone.popSize;
            clone=(T)clone.next;
            iClone++;
        }
        CloneCounts.add(ret);
    }

//    /**
//     * writes out the parent ids to a file
//     * @param parentIDsOut the file to write to
//     * @param delim delimeter to separate each id entry
//     */
//    public void WriteParentIDs(FileIO parentIDsOut,String delim){
//        for (int id : parentIDs) {
//            if(parentIDsOut.binary){
//                parentIDsOut.WriteBinInt(id);
//            }
//            else{
//                parentIDsOut.Write(id+delim);
//            }
//        }
//    }
//
//    /**
//     * writes the mutation info strings to a file
//     * @param mutationInfoOut the file to write to
//     * @param delim delimeter to separate each info entry
//     */
//    public void WriteMutationInfo(FileIO mutationInfoOut,String delim){
//        for(String info: allGenomeInfos){
//            if(mutationInfoOut.binary){
//                mutationInfoOut.WriteBinString(info);
//            }
//            else{
//                mutationInfoOut.Write(mutationInfoOut+delim);
//            }
//        }
//    }

    /**
     * writes the mutation info of all clone types that are currently alive in the form id delim info delim id delim info delim... to a file
     * @param mutationInfoOut file to write to
     * @param delim delimiter to separate each entry
     */
    public void WriteMutationInfoLiving(FileIO mutationInfoOut, String delim){
        T curr=livingCloneInfos;
        while(curr!=null) {
            mutationInfoOut.Write(curr.id+delim+curr.GenomeInfoStr()+delim);
            curr=(T)curr.next;
        }
    }

    /**
     * gets a concatenated string of the full lineage of the particular clone in the form info delim info delim...
     * @param id id of genome of interest
     * @param delim used to separate genome strings along lineage
     */
    public String FullLineageInfoStr(int id, String delim){
        StringBuilder sb=new StringBuilder();
        while(id!=-1){
            sb.append(allGenomeInfos.get(id)+delim);
            id=parentIDs.get(id);
        }
        return sb.toString();
    }

    /**
     * writes the full lineage of all clone types that are currently alive in the form id outerdelim info innerdelim info innerdelim... outerdelim id outerdelim info...
     * @param lineageInfoOut the file to write to
     * @param innerDelim used to separate info entries
     * @param outerDelim used to separate ids from info
     */
//    public void WriteAllLineageInfoLiving(FileIO lineageInfoOut,String innerDelim,String outerDelim){
//        T curr=livingCloneInfos;
//        while(curr!=null){
//            if(lineageInfoOut.binary){
//                lineageInfoOut.WriteBinString(curr.id+innerDelim+ FullLineageInfoStr(curr.id,innerDelim)+outerDelim);
//            }
//            else {
//                lineageInfoOut.Write(curr.id+innerDelim+ FullLineageInfoStr(curr.id, innerDelim) + outerDelim);
//            }
//            curr=(T)curr.next;
//        }
//    }

//    /**
//     * writes the full lineage of all clone types that have ever existed in the form outerdelim info innerdelim info innerdelim... outerdelim info...
//     * @param lineageInfoOut the file to write to
//     * @param innerDelim used to separate info entries
//     * @param outerDelim used to separate lineages
//     */
//    public void WriteAllLineageInfo(FileIO lineageInfoOut,String innerDelim,String outerDelim) {
//        for (int i = 0; i < parentIDs.size(); i++) {
//            if(lineageInfoOut.binary){
//                lineageInfoOut.WriteBinString(FullLineageInfoStr(i,innerDelim)+outerDelim);
//            }
//            else{
//                lineageInfoOut.Write(i + innerDelim + FullLineageInfoStr(i,innerDelim)+outerDelim);
//            }
//        }
//    }

    /**
     * writes all clone population counts stored to a file in the form outerdelim id innerdelim info innerdelim id innerdelim info... outerdelim...
     * @param cloneCountsOut file to write clonecounts to
     * @param innerDelim used to separate ids and infos
     * @param outerDelim used to separate clonecount events
     */
    public void WriteClonePops(FileIO cloneCountsOut, String innerDelim, String outerDelim){
        for (int[] counts : CloneCounts) {
            cloneCountsOut.WriteDelimit(counts,innerDelim);
            cloneCountsOut.Write(outerDelim);
        }
    }
}