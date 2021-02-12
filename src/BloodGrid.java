import Epigenome.EpigenomeSkeleton;
import Epigenome.GenomeTracker;
import HAL.GridsAndAgents.AgentGrid0D;
import HAL.Rand;
import HAL.Tools.FileIO;

import java.text.DecimalFormat;
import java.util.ArrayList;

public class BloodGrid extends AgentGrid0D<BloodCell> {
    GenomeTracker<Epigenome> EpiGenomeStorage;
    Rand RN=new Rand();
    DecimalFormat df = new DecimalFormat("0.0000");

    double DiseaseProp;

    /*
    Initialization Variables
     */
    final double[] data=new double[]{0.000143149,0.000269805,0.00091566,0.00114652,0.000772511,0.000477827,0.000380739,0.000399495,0.000492114,0.000644973,0.000851174,0.00116682,0.001669415,0.002409655,0.003492571,0.005188465,0.007823022,0.011828891,0.018141149,0.027322078,0.039480516,0.053155456,0.065228666,0.074181935,0.080673419,0.084714151,0.086727677,0.085985285,0.081905829,0.073437172,0.059573462,0.043273949,0.028662001,0.01820751,0.011455875,0.00739236,0.004952089,0.003490475,0.002512922,0.001888525,0.001472206,0.001215529,0.001038454,0.000912295,0.000827674,0.000718285,0.000575135,0.00039293,0.000226612,0.000140612,4.50E-05};
//    final double[] data=GetData(); // TODO use this to get from File
    int[][] BinnedEpiSites;
    double[] ecdfProbs = new double[]{0.0196078431372549,0.0392156862745098,0.0588235294117647,0.0784313725490196,0.0980392156862745,0.117647058823529,0.137254901960784,0.156862745098039,0.176470588235294,0.196078431372549,0.215686274509804,0.235294117647059,0.254901960784314,0.274509803921569,0.294117647058824,0.313725490196078,0.333333333333333,0.352941176470588,0.372549019607843,0.392156862745098,0.411764705882353,0.431372549019608,0.450980392156863,0.470588235294118,0.490196078431373,0.509803921568627,0.529411764705882,0.549019607843137,0.568627450980392,0.588235294117647,0.607843137254902,0.627450980392157,0.647058823529412,0.666666666666667,0.686274509803922,0.705882352941177,0.725490196078431,0.745098039215686,0.764705882352941,0.784313725490196,0.803921568627451,0.823529411764706,0.843137254901961,0.862745098039216,0.882352941176471,0.901960784313726,0.92156862745098,0.941176470588235,0.96078431372549,0.980392156862745,1};;
//    double[] ecdfProbs = GetDataECDF(); // TODO use this to get from File
    double[][] errorProbs = new double[][]{{0,2.00E-05,0.00004,0.00006,0.00008,1.00E-04,0.00012,0.00014,0.00016,0.00018,0.0002,0.00022,0.00024,0.00026,0.00028,0.0003,0.00032,0.00034,0.00036,0.00038,0.0004,0.00042,0.00044,0.00046,0.00048,0.0005,0.00052,0.00054,0.00056,0.00058,0.0006,0.00062,0.00064,0.00066,0.00068,0.0007,0.00072,0.00074,0.00076,0.00078,0.0008,0.00082,0.00084,0.00086,0.00088,0.0009,0.00092,0.00094,0.00096,0.00098,0.001},{0.001,0.00098,0.00096,0.00094,0.00092,0.0009,0.00088,0.00086,0.00084,0.00082,0.0008,0.00078,0.00076,0.00074,0.00072,0.0007,0.00068,0.00066,0.00064,0.00062,0.0006,0.00058,0.00056,0.00054,0.00052,0.0005,0.00048,0.00046,0.00044,0.00042,0.0004,0.00038,0.00036,0.00034,0.00032,0.0003,0.00028,0.00026,0.00024,0.00022,0.0002,0.00018,0.00016,0.00014,0.00012,0.0001,0.00008,6.00E-05,0.00004,0.00002,0}};


    /*
    Percent Divisions
     */
    ArrayList<Integer> Divisions=new ArrayList<>();
    int divs=0;
    ArrayList<Integer> Losses=new ArrayList<>();
    int loss=0;

    /*
    Disease information
     */
    boolean DiseaseInit=false;

    public BloodGrid(EpigenomeSkeleton epiArray){
        super(BloodCell.class);
        BinnedEpiSites = AssignEpiSiteBins();

        EpiGenomeStorage = new GenomeTracker<>(new Epigenome(1f,0f,1, epiArray.InitializeMethylationState(RN), errorProbs, RN, BinnedEpiSites), false, false);

        Init(epiArray);
        System.out.println("Initialization Completed");
    }

    public void RunStep(){
        divs=0;
//        System.out.println("Initial Population Size: " + GetPop());
//        System.out.println("Initial Genome Mutants: " + EpiGenomeStorage.NumActiveMutantTypes());
        if(CON.DISEASE == true){
            SetDiseasePop();
        }

        for(BloodCell c: this) {
            if(CON.DISEASE==true & CON.DISEASEINIT==GetTick() & DiseaseInit==false){
                c.CellStatus=1;
                DiseaseInit=true;
            }
            CheckDivision(c);
//            if(c.DivDaysLeft<=0){
//                if(CON.DISEASE == true && c.CellStatus==0){
//                    loss+=1; // Symmetric Loss
//                    c.myEpigenome.epiData = null;
//                    c.myEpigenome.DisposeClone();
//                    c.Dispose();
//                }
//            }
        }

        Losses.add(loss);
        Divisions.add(divs);
        CleanShuffle(RN);
        IncTick();
    }

    public void SetDiseasePop(){
        double val=0;
        for(BloodCell c: this) {
            if(c.CellStatus==1){
                val+=1;
            }
        }
        DiseaseProp = val/GetPop();
    }

    public void CheckDivision(BloodCell c){
        if(CON.DISEASE == true & c.CellStatus==1 & DiseaseProp<=CON.CARRYCAP & RN.rn.Double()<CON.DISEASEDIVPROB){ // Disease state
                BloodCell d1 = NewAgent(); // Symmetric division of stem cells in disease...No Loss.
                d1.BirthClonalPop(c.cloneID, c.myEpigenome.NewChild().PossiblyMutate(), CON.DAYS*2, c.CellStatus);
                divs+=1;
        } else if(CON.DISEASE == true & c.CellStatus==1 & DiseaseProp>CON.CARRYCAP) {
            NormalDivCase(c);
        } else {
            NormalDivCase(c);
        }
    }

    public void NormalDivCase(BloodCell c) {
        if (c.DivDaysLeft > 0 & RN.rn.Double() < CON.DIVPROP) {
            BloodCell d = NewAgent(); // One Daughter is born, the other is lost (thus only one is kept as stem, the new one)
            d.BirthClonalPop(c.cloneID, c.myEpigenome.NewChild().PossiblyMutate(), c.DivDaysLeft, c.CellStatus);
            if (c.myEpigenome.mutated) {
                c.myEpigenome.epiData = null;
                c.myEpigenome.DisposeClone();
            }
            c.Dispose();
            divs += 1;
            loss += 1; // Assymetric division loss
        }
    }

    private void Init(EpigenomeSkeleton epiArray){
        BloodCell c = NewAgent();
        c.BirthClonalPop(100, EpiGenomeStorage.NewProgenitor(), GetDivDays(), 0);
        for (int i = 0; i < CON.INITPOP; i++) {
            BloodCell d = NewAgent();
            d.BirthClonalPop(i, c.myEpigenome.NewChild().PossiblyMutate(), GetDivDays(), 0);
            SetInitialEpiGenome(d);
        }
        if(c.myEpigenome.mutated){
            c.myEpigenome.DisposeClone();
        }
        c.Dispose();
        CleanShuffle(RN);
    }

    private int GetDivDays(){
        return(RN.rn.Int(CON.STEMDIVSALLOWED[1]-CON.STEMDIVSALLOWED[0])+CON.STEMDIVSALLOWED[0]);
    }

    public void ReportPopulation(){
        System.out.println(GetTick() + "\t" + GetPop());
    }

    private void SetInitialEpiGenome(BloodCell c){
        int methylVal;
        int binid;
        for (int i = 0; i < CON.EPISITES; i++) {
            binid = BinnedEpiSites[0][i];
            methylVal = RN.Binomial(1,ecdfProbs[binid]);
            c.myEpigenome.epiData.Set(0,i, methylVal==1 ? true:false);
            methylVal = RN.Binomial(1,ecdfProbs[binid]);
            c.myEpigenome.epiData.Set(1,i, methylVal==1 ? true:false);
        }
    }

    private double[] GetData(){
        ArrayList<Double> vals = ParseData();

        double[] val = new double[vals.size()];
        for (int i = 0; i < vals.size(); i++) {
            val[i] = vals.get(i);
        }
        return(val);
    }

    private ArrayList<Double> ParseData(){
        FileIO reader = new FileIO(System.getProperty("user.dir") + "/src/data/EpiAlleleTable.23Apr20.txt", "r");
        String[] line = reader.ReadLineDelimit("\t");
        line = reader.ReadLineDelimit("\t");
        ArrayList<Double> percentages = new ArrayList<>();

        while(line != null){
            percentages.add(Double.parseDouble(line[0]));
            line = reader.ReadLineDelimit("\t");
        }
        return(percentages);
    }

//    private ArrayList<Double> ParseRate(){
//        FileIO reader = new FileIO(System.getProperty("user.dir") + "/src/data/EpiAlleleTable.23Apr20.txt", "r");
//        String[] line = reader.ReadLineDelimit("\t");
//        line = reader.ReadLineDelimit("\t");
//        ArrayList<Double> rate = new ArrayList<>();
//
//        while(line != null){
//            rate.add((Double.parseDouble(line[4]) + Double.parseDouble(line[5])) / 2f);
//            line = reader.ReadLineDelimit("\t");
//        }
//        return(rate);
//    }

    private double[] GetDataECDF(){
        ArrayList<Double> vals = ParseDataECDF();
        double[] val = new double[vals.size()];
        for (int i = 0; i < vals.size(); i++) {
            val[i] = vals.get(i);
        }
        return(val);
    }

    private ArrayList<Double> ParseDataECDF(){
        FileIO reader = new FileIO(System.getProperty("user.dir") + "/src/data/ECDF.ready.txt", "r");
        String[] line = reader.ReadLineDelimit("\t");
        line = reader.ReadLineDelimit("\t");
        ArrayList<Double> percentages = new ArrayList<>();

        while(line != null){
            percentages.add(Double.parseDouble(line[0]));
            line = reader.ReadLineDelimit("\t");
        }
        return(percentages);
    }

    public double[] MethylationResults(){
        double[] methylationLevel = new double[CON.EPISITES];
        double[] ones = new double[CON.EPISITES];

        for(BloodCell c: this){
            for (int i = 0; i < CON.CHROMNUM*2; i++) {
                for (int j = 0; j < CON.EPISITES; j++) {
                    ones[j]+=c.myEpigenome.epiData.GetInt(i,j);
                }
            }
        }

        for (int i = 0; i < ones.length; i++) {
            methylationLevel[i] = ones[i]/(double)(GetPop()*(CON.CHROMNUM*2));
        }

        return(methylationLevel);
    }

    public int[][] AssignEpiSiteBins(){
        int[][] finalRet = new int[2][];
        int[] ret = new int[CON.EPISITES];
        int[] binned = new int[data.length];

        RN.Multinomial(data, CON.EPISITES, binned); // Sets values within binned with number of episites in that bin

        int idx = 0;
        for (int i = 0; i < binned.length; i++) {
            for (int j = 0; j < binned[i]; j++) {
                ret[idx] = i;
                idx+=1;
            }
        }
        finalRet[0] = ret;
        finalRet[1] = binned;

        return(finalRet);
    }

    public double[][] MethylationResults2(){
        double[][] ones = new double[CON.INITPOP][CON.EPISITES];

        int cell = 0;
        for(BloodCell c: this){
            for (int i = 0; i < CON.CHROMNUM*2; i++) {
                for (int j = 0; j < CON.EPISITES; j++) {
                    ones[cell][j]+=c.myEpigenome.epiData.GetInt(i,j);
                }
            }
            cell +=1;
        }

        return(ones);
    }

    public ArrayList MethylationWritableResults2(){
        double[][] methylationResults = MethylationResults2();

        ArrayList<String> methylationWritable = new ArrayList<>();

        String header = "MethylSite\t";
        for (int i = 0; i < CON.INITPOP; i++) {
            header+="c"+i+"\t";
        }
        header+="\n";
        methylationWritable.add(header);

        for (int j = 0; j < CON.EPISITES; j++) {
            String bodyline = Integer.toString(j);
            for (int i = 0; i < CON.INITPOP; i++) {
                bodyline += "\t" + methylationResults[i][j];
            }
            bodyline+="\n";
            methylationWritable.add(bodyline);
        }
        return(methylationWritable);
    }

    public ArrayList MethylationWritableResults(){
        double[] methylationResults = MethylationResults();

        ArrayList<String> methylationWritable = new ArrayList<>();

        methylationWritable.add("MethylSite\tMethylationLevel\n");
        for (int i = 0; i < methylationResults.length; i++) {
            methylationWritable.add(i + "\t" + methylationResults[i] + "\n");
        }
        return(methylationWritable);
    }

    public ArrayList GetBinnedResults(){
        ArrayList<String> ret = new ArrayList<>();

        ret.add("Bin\tCpGLocis\n");
        for (int i = 0; i < BinnedEpiSites[1].length; i++) {
            ret.add(i + "\t" + BinnedEpiSites[1][i] + "\n");
        }

        return(ret);
    }

    public ArrayList GetECDFResults(){
        ArrayList<String> ret = new ArrayList<>();

        ret.add("Bin\tECDF\n");
        for (int i = 0; i < ecdfProbs.length; i++) {
            ret.add(i + "\t" + ecdfProbs[i] + "\n");
        }

        return(ret);
    }

    public void ReportDivs(){
        int val = 0;
        for (int i = 0; i < Divisions.size(); i++) {
            val += Divisions.get(i);
        }
        System.out.println("Daily Turnover (%): " + ((val/(double)Divisions.size())/CON.INITPOP) * 100. );
    }

    public void ReportLosses(){
        int val = 0;
        for (int i = 0; i < Losses.size(); i++) {
            val += Losses.get(i);
        }
        System.out.println("Daily Losses (%): " + ((val/(double)Divisions.size())/CON.INITPOP) * 100. );
    }

    public String SamplePop(){
        int diseasePop=0;
        int healthyPop=0;
        for(BloodCell c : this){
            if(c!=null){
                if(c.CellStatus==0){
                    healthyPop+=1;
                } else if (c.CellStatus==1){
                    diseasePop+=1;
                }
            }
        }
        return((GetTick() + "\t" + GetPop() + "\t" + df.format(healthyPop/(double)GetPop()) + "\t" + df.format(diseasePop/(double)GetPop())));
    }

    public double[] GetAbNormalPopProp(){
        int diseasePop=0;
        int healthyPop=0;
        for(BloodCell c : this){
            if(c!=null){
                if(c.CellStatus==0){
                    healthyPop+=1;
                } else if (c.CellStatus==1){
                    diseasePop+=1;
                }
            }
        }
        return(new double[]{healthyPop/(double)GetPop(),diseasePop/(double)GetPop()});
    }

    public double GetVariance(){
        double[] d = MethylationResults();
        double num=0;
        double xi=mean(d);
        for (int i = 0; i < d.length; i++) {
            num+=Math.pow((d[i]-xi),2);
        }
        return(num/(d.length-1));
    }

    public double mean(double[] d){
        double sum=0;
        for (int i = 0; i < d.length; i++) {
            sum+=d[i];
        }
        return(sum/d.length);
    }

}

