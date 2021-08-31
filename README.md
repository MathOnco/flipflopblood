## Summary

Simulations of hematopoietic stem cells to accompany our paper "Cell lineage tracing with molecular clocks based on fluctuating DNA methylation."

The provided ```.jar``` file is used to exectue the model while the two accompanying files (```run_simulations.sh``` and ```GeneratePlots.py```) allow you to recreate the plots used within the paper.

## Dependencies

1. Python (>=3.6)
2. [Java 1.8](https://www.java.com/en/download/manual.jsp)
3. [Java JDK](https://www.oracle.com/java/technologies/javase-downloads.html)

## Execution

1. Clone the repository: ```git clone https://github.com/MathOnco/ticktockblood.git```

2. Navigate into the directory from the commandline ```cd ticktockblood```

3. Execute the simulations ```bash run_simulations.sh```

4. Once simulations have completed, run the python script to generate figures using the following command ```python GeneratePlots.py Normal 2500 ALL.0.9 300 ALL.0.5 300 ALL.0.2 300 CL.0.9 1500 CL.0.5 1500 CL.0.2 1500 CHIP.0.9 2500 CHIP.0.5 2500 CHIP.0.2 2500```

This will produce two ```.png``` files from the simulations.

## Additional informaiton

To execute the model using any parameter set you want refer to the help menu from the jar file:
```
java -jar ./TickTockBlood.jar --help
```
Source code is located within the ```./src/``` directory. This includes dependencies within the HAL framework and those required for the jar file build.
