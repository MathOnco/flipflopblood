import HAL.GridsAndAgents.Agent0D;

public class BloodCell extends Agent0D<BloodGrid> {
    Epigenome myEpigenome;
    int cloneID;
    int StemStatus=1;
    int DivDaysLeft;
    int CellStatus;

    public void BirthClonalPop(int id, Epigenome epigenome, int DivDays, int CellDiseaseStatus){
        cloneID=id;
        myEpigenome=epigenome;
        DivDaysLeft=DivDays;
        CellStatus=CellDiseaseStatus;
    }

}
