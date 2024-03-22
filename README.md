# DoublyValidDoublySharpSensitivity

Replication code for "Doubly-Valid/Doubly-Sharp Sensitivity Analysis for Causal Inference with Unmeasured Confounding" by Jacob Dorn, Kevin Guo, and Nathan Kallus (https://arxiv.org/abs/2112.11449)

There are two separate tasks: one to replicate the results in Section 7.1 and one for Section 7.2. The two projects share the file `DVDS.R` and `simplified_dvds.R`, but otherwise can be run independently.

The file `DVDS.R` is a set of underlying functions used in estimation. The file `simplified_dvds.R` are wrapper DVDS and ZSB estimation functions for binary and continuous outcomes. 

The first task is the simulations (Section 7.1). To reproduce the figures in the paper: 
1. First run `dvds_sim.R`.  This may take one or several hours, and will output data in the folder `/Simulations/dvds_output/`
2. Then run `dvds_plotting.R`.  This is roughly instantaneous and will output plots in `/Simulations/dvds_output/`

The second task is to analyze the RHC data (Section 7.2). The file `run_DVDS_RHCExploration.R` runs the real data application using the data in `/Real_data/rhc.csv`. The script generates outputs in `/Real_data/Results/`. Another file called temp is generated but not substantive. 

The results reported in the paper were run on an Apple MacBook Pro running R version 4.2.1. Last updated: March 21, 2024.