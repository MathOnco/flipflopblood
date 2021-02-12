import sys
import matplotlib.pyplot as plt
from scipy.stats import sem, t
from scipy import mean
import numpy as np

def ci(tmpt_arr, conf=0.95):
    m = mean(tmpt_arr)
    se = sem(tmpt_arr)
    h = se * t.ppf((1 + conf) / 2, len(tmpt_arr) - 1)
    return((m-h,m+h))

def getBetaDistPlots(f):
    props = ["Normal_1","ALL.0.9_1","CL.0.9_1","CHIP.0.9_1"]

    fig, axs = plt.subplots(1, len(props), sharey=False, tight_layout=True, figsize=(8, 3), dpi=300)

    for i, sam in enumerate(props):
        with open("./Output.betadist.%s.csv"%(props[i]), "r") as inFile:
            lines = inFile.readlines()
        del lines[0]
        lines = [item.replace("\n","").split("\t") for item in lines]
        lines = [float(item[len(item)-1]) for item in lines]

        axs[i].hist(lines, density=True, bins=71)
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['top'].set_visible(False)
        axs[i].set_xlabel("Beta")
        axs[i].set_ylabel("Probability Density")

    plt.savefig('betaDists.png', bbox_inches='tight')

def getVariancePlot(f, tmpts, reps=6):
    linetypes = ['dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    color = ['#000000','#0061a6','#f97523','#05853c','#0061a6','#f97523','#05853c','#0061a6','#f97523','#05853c']

    fig, ax = plt.subplots()
    for i, sam in enumerate(f):
        replicate_data = np.zeros((int(tmpts[i])+1, reps-1))
        for j in range(1,reps):
            with open("./Output.variance.%s_%s.csv"%(sam,j), "r") as inFile:
                lines = inFile.readlines()
            divs = [int(item.replace("\n","").split(",")[0]) for item in lines]
            vals = [float(item.replace("\n","").split(",")[1]) for item in lines]

            for k in range(0,tmpts[i]+1):
                replicate_data[k][j-1]=vals[k]

        # get mean and ci at each timepoint
        sim_mean = []
        sim_ci_low = []
        sim_ci_high = []
        for k in range(0,tmpts[i]+1):
            ci_tmpt = ci(replicate_data[k])
            sim_mean.append(mean(replicate_data[k]))
            sim_ci_low.append(ci_tmpt[0])
            sim_ci_high.append(ci_tmpt[1])

        ax.plot(divs,sim_mean, linestyle=linetypes[i], color=color[i], linewidth=1, label=sam)
        ax.fill_between(divs, sim_ci_high, sim_ci_low, facecolor=color[i], alpha=0.5)
        ax.set(xlim=(0, 2500))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    leg = ax.legend()
    plt.xlabel("Simulated Divisions from Initiation")
    plt.ylabel("Variance of the Beta Distribution\nof Oscillatory CpG Loci")
    plt.savefig('variances.png', bbox_inches='tight')

def main():
    args = sys.argv
    del args[0]

    f = args[::2]
    tmpts = args[1::2]
    tmpts = [int(val) for val in tmpts]

    getVariancePlot(f, tmpts)
    getBetaDistPlots(f)

if __name__=="__main__":
    main()
