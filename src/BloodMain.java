import HAL.Tools.FileIO;
import picocli.CommandLine;
import Epigenome.EpigenomeSkeleton;
import Epigenome.IntArrayEpigenome;

import java.util.ArrayList;

class CON{
    /*
    Global Parameters
     */
    static int DAYS=500; // Simulation time. 1 day == 1 division for each normal HSC
    static int INITPOP=20; // initialization population
    static int[] STEMDIVSALLOWED=new int[]{101100,101101}; // min, max. Number of divisions a single HSC is permitted before ordered replacement
    static double DIVPROP=1; // 0.5 provides the variance level appropriate for 0.015 observed in Normal.

    /*
    Epigenome Parameters
     */
    static int CHROMNUM=1; // Number of chromosomes (paired means this number times two) # No need to adjust this
    static int EPISITES=10000; // Number of sites on each chromosome. Proper Number 27634

    /*
    Disease States
     */
    static boolean DISEASE=true; // Set to true to run disease simulations
    static double DISEASEDIVPROB=0.1; // ALL=0.1, Chronic Leukemia=0.005, CHIP=0.00225?
    static int DISEASEINIT=10;
    static double CARRYCAP=0.2;

    static boolean JAR=true;

    static String RUNNAME="ALL.0.2_5";
}

@CommandLine.Command(name = "Ticktock Blood Model",
        header="Model execution can be carried out using several specified options allowing for rapid deployment options.",
        mixinStandardHelpOptions = false,
        showDefaultValues = true,
        sortOptions = false,
        headerHeading = "@|bold,underline,red Usage|@:%n%n",
        synopsisHeading = "%n",
        descriptionHeading = "%n@|bold,underline,red Description|@:%n%n",
        parameterListHeading = "%n@|bold,underline,red Parameters|@:%n",
        optionListHeading = "%n@|bold,underline,red Options|@:%n",
        version = "Blood",
        description = "")
public class BloodMain implements Runnable {
    @CommandLine.Option(names = { "-V", "--version" }, description = "Display version and exit.") boolean versionRequested;
    @CommandLine.Option(names = { "-h", "--help" }, usageHelp = true, description = "Display a help message and exit.") boolean helpRequested;
    @CommandLine.Option(names = { "-t", "--runtime"}, description="Time to execute model for.") int days=CON.DAYS;
    @CommandLine.Option(names = { "-s", "--popsize"}, description="Initial Population Size.") int pop=CON.INITPOP;
    @CommandLine.Option(names = { "-d", "--divprop"}, description="Division Probability of Normal") double divprop=CON.DIVPROP;
    @CommandLine.Option(names={"-chr","--chroms"}, description="Number of chromosomes (x2).") int chroms=CON.CHROMNUM;
    @CommandLine.Option(names={"-chrs","--methylationSites"}, description="Number of methylation sites.") int methSites=CON.EPISITES;
    @CommandLine.Option(names = { "-a", "--disease"}, description="Disease Model. Boolean.") boolean disease=CON.DISEASE;
    @CommandLine.Option(names = { "-c", "--expansionrate"}, description="Division Probability of Cancer") double expansionrate=CON.DISEASEDIVPROB;
    @CommandLine.Option(names = { "-i", "--inductionday"}, description="Inducition day for disease.") int diseaseday=CON.DISEASEINIT;
    @CommandLine.Option(names = { "-p", "--carry"}, description="Carrying cap 0-1.") double carryingcap=CON.CARRYCAP;
    @CommandLine.Option(names={"-rn","--runname"}, description="Name for the simulation outputs.") String runName=CON.RUNNAME;

    public void run(){
        assert !helpRequested;
        CON.DAYS=days;
        CON.INITPOP=pop;
        CON.DIVPROP=divprop;
        CON.CHROMNUM=chroms;
        CON.EPISITES=methSites;
        CON.DISEASE=disease;
        CON.DISEASEDIVPROB=expansionrate;
        CON.DISEASEINIT=diseaseday;
        CON.CARRYCAP=carryingcap;
        CON.RUNNAME=runName;

        System.out.println(CON.DAYS);

        System.out.println(runName);
        RunModel();
    }

    static double[][] methylationOut;
//    static ArrayList<String> RunStatsOut = new ArrayList<>();
    static EpigenomeSkeleton epiArray;
    static BloodGrid G;
    static double[] varianceArray;

    public static void main(String[] args) {
            if(CON.JAR){
                CommandLine.run(new BloodMain(), args);
            } else {
                RunModel();
            }
    }

    public static void RunModel() {
        methylationOut = new double[CON.DAYS/200][CON.EPISITES];
        epiArray = new IntArrayEpigenome(CON.CHROMNUM*2,CON.EPISITES);
        varianceArray = new double[CON.DAYS+1];

        System.out.println("IDE Run");
        G = new BloodGrid(epiArray);

        int count=0;
        RecordMethylation(G.GetTick(), G.MethylationResults());
        varianceArray[G.GetTick()]=G.GetVariance();

        System.out.println("Beginning Runsteps.");
        while(G.GetTick()<CON.DAYS){

            G.RunStep();

            varianceArray[G.GetTick()]=G.GetVariance();

            if(G.GetTick()%200==0){
                String Stats=G.SamplePop();
//                RunStatsOut.add(Stats);
                System.out.println(Stats);
                RecordMethylation(count, G.MethylationResults());
                count +=1;
                System.out.println("Memory (Gb): " + Runtime.getRuntime().totalMemory()/1073741824.0);
            }
        }
        System.out.println("Memory (Gb): " + Runtime.getRuntime().totalMemory()/1073741824.0);
        System.out.println("Complete. Writing outputs...");
        writeResultsGeneric(varianceArray, "./Output.variance." + CON.RUNNAME + ".csv");
//        writeRunStats(RunStatsOut, "Output.RunStats.txt");
        writeResults(PrepareMethylationOutput(), "./Output.betadist." + CON.RUNNAME + ".csv");
//        writeResults(G.GetBinnedResults(), "Output.Bins.txt");
//        writeResults(G.GetECDFResults(), "Output.ECDF.txt");

    }

    public static void writeResultsGeneric(double[] arrIn, String theFileName){
        FileIO fileWriter = new FileIO(theFileName, "w");
        for (int i = 0; i < arrIn.length; i++) {
            fileWriter.Write(Integer.toString(i) + "," + Double.toString(arrIn[i]) + "\n");
        }
        fileWriter.Close();
    }

    public static void writeResults(ArrayList<String> methylationOutput, String theFileName){
        FileIO mylationWriter = new FileIO(theFileName, "w");
        for (int i = 0; i < methylationOutput.size(); i++) {
            mylationWriter.Write(methylationOutput.get(i));
        }
        mylationWriter.Close();
    }
//
//    public static void writeRunStats(ArrayList<String> RunStats, String theFileName){
//        FileIO mylationWriter = new FileIO(theFileName, "w");
//        for (int i = 0; i < RunStats.size(); i++) {
//            mylationWriter.Write(RunStats.get(i) + "\n");
//        }
//        mylationWriter.Close();
//    }

    public static void RecordMethylation(int tick, double[] meth){
        for (int i = 0; i < meth.length; i++) {
            methylationOut[tick][i] = meth[i];
        }
    }

    public static ArrayList PrepareMethylationOutput(){
        ArrayList<String> methylationWritable = new ArrayList<>();

        String header="MethylSite\t";
        for (int i = 0; i < methylationOut.length; i++) {
            if(i<methylationOut.length-1){
                header+=i + "\t";
            } else {
                header+=i + "\n";
            }
        }
        methylationWritable.add(header);
        for (int i = 0; i < methylationOut[0].length; i++) {
            String site=i + "\t";
            for (int j = 0; j < methylationOut.length; j++) {
                if(j!=methylationOut.length-1){
                    site+=methylationOut[j][i] + "\t";
                } else {
                    site+=methylationOut[j][i] + "\n";
                }
            }
            methylationWritable.add(site);
        }

        return(methylationWritable);
    }

}

