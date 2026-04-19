# A Rule-Based Modeling Approach for Studying AMPA Receptors at the Synapse

## Introduction of the repository

The original and modified Kappa code, results, and visualization files are available here. 

**ATTENTION**: The results related to **04_wLTP_burst** are provided in the .zip file uploaded to Blackboard Learn as they are *too large* to upload to github.com. These three .csv files should be copied to the **“Result”** folder for successful visualization. The usage of the Kappa code is explained in the “Instruction_for_running.txt” file in the **“Code”** folder. 

## Overview

In this project, a rule-based modeling language *Kappa* was used to simulate molecular interactions in three pre-constructed AMPAR models above using five LTP induction protocols, including single theta burst train (siTBS), compressed theta burst (cTBS), spaced theta burst (spTBS), weak LTP burst (wLTP), and high frequency stimulation (HFS) *(Bliss and Lømo, 1973; Grover et al., 2009; Park et al., 2016)*. We aim to characterize the dynamics of AMPAR trafficking and abundance under different LTP stimulation conditions *in silico*, guiding future wet-lab experiments.

## Method summary

The Kappa code simulating the dynamics of AMPARs at the synapse in the Slot, Receptor, and Exocytosis models was adapted from the original models provided by the instructor. In this project, we mainly focused on adding five different LTP induction protocols mainly by modifying the **“Modifications.ka” files** in each model to implement the different LTP induction protocols. **The modified code is clearly documented at the end of the "Modifications.ka" files.** For each model under each LTP protocol, 0s to 1s was set as the baseline, and calcium stimulation started at 1.0s per the protocol. After the final stimulation, we observed the effect of the LTP protocol on the variables for at least 3 seconds. The detailed time points for each LTP protocol, during and after stimulation, are shown in Supplementary Fig. 1B-C in *Supporting materials*.

After modifying the code, the models were run from the command line using KaSim, as described in the user manual at https://kappalanguage.org/. For each simulation, KaSim was provided with the model files in sequence using repeated -i inputs: Agents.ka, Variables.ka, 8020.ka, Rules.ka, Observables.ka, and Modifications.ka. For siTBS, cTBS, wLTP and HFS, outputs were recorded every 0.0001 s; for spTBS, outputs were recorded every 0.005 s because of the much longer simulation window. The total simulation time was set according to the LTP protocol shown in Supplementary Fig. 1C in *Supporting materials*. For details please look up the “Instruction_for_running.txt” file in the **“Code”** folder. 

## Repository structure and file explanation
```
├─Code
│  ├─00_Original_Model
│  │  ├─Exocytosis
│  │  ├─Receptor
│  │  └─Slot
│  ├─01_Single_theta_burst_train
│  │  ├─Exocytosis
│  │  ├─Receptor
│  │  └─Slot 
│  ├─02_Compressed_theta_burst
│  │  ├─Exocytosis
│  │  ├─Receptor
│  │  └─Slot
│  ├─03_Spaced_theta_burst
│  │  ├─Exocytosis
│  │  ├─Receptor
│  │  └─Slot
│  ├─04_wLTP_burst
│  │  ├─Exocytosis
│  │  ├─Receptor
│  │  └─Slot
│  └─05_HFS
│      ├─Exocytosis
│      ├─Receptor
│      └─Slot
├─Result
└─Visualisation
```

### 1. Code/

Contains the original and modified Kappa code. The usage of the Kappa code is explained in the “Instruction_for_running.txt” file in the **“Code”** folder. 

- **00_Original_Model**: The original Kappa code simulating the movement and trafficking of AMPARs in the posy-synaptic cells given by the instructor.

- **01_Single_theta_burst_train**: The Kappa code simulating the movement and trafficking of AMPARs under single theta burst train stimulation in three models.

- **02_Compressed_theta_burst**: The Kappa code simulating the movement and trafficking of AMPARs under compressed theta burst train stimulation in three models.

- **03_Spaced_theta_burst**: The Kappa code simulating the movement and trafficking of AMPARs under single theta burst train stimulation in three models.

- **04_wLTP_burst**: The Kappa code simulating the movement and trafficking of AMPARs under weak LTP burst stimulation in three models.

- **05_HFS**: The Kappa code simulating the movement and trafficking of AMPARs under high frequency stimulation in three models.


### 2. Result/

Contains 12 CSV files (the other 3 CSV files related to the results of wLTP burst were uploaded to Blackboard Learn), which are results after the models were run from the command line using KaSim.

### 3. Visualisation/

Contains two R files (Figure 1.R & Figure 2.R) that could visualize the data in the Result folder to generate Figure 1&2 in the main text. Before running the script, please make sure the file path is correctly defined and relevant R packages are installed with the correct version.

R session info:
```
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] patchwork_1.2.0 purrr_1.0.4     tidyr_1.3.1     readr_2.1.5     stringr_1.5.1  
[6] dplyr_1.1.4     ggplot2_4.0.2  

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5        cli_3.6.3          rlang_1.1.4        stringi_1.8.4     
 [5] generics_0.1.3     S7_0.2.1           glue_1.7.0         hms_1.1.3         
 [9] scales_1.4.0       fansi_1.0.6        grid_4.4.1         tibble_3.2.1      
[13] tzdb_0.5.0         lifecycle_1.0.4    compiler_4.4.1     RColorBrewer_1.1-3
[17] pkgconfig_2.0.3    rstudioapi_0.16.0  farver_2.1.2       R6_2.5.1          
[21] tidyselect_1.2.1   utf8_1.2.4         pillar_1.9.0       magrittr_2.0.3    
[25] tools_4.4.1        withr_3.0.1        gtable_0.3.6 
```

## Reference

Bliss, T.V.P. & Lømo, T. (1973) ‘Long-lasting potentiation of synaptic transmission in the dentate area of the anaesthetized rabbit following stimulation of the perforant path’, The Journal of Physiology, 232, 331–356.

Grover, L.M., Kim, E., Cooke, J.D. & Holmes, W.R. (2009) ‘LTP in hippocampal area CA1 is induced by burst stimulation over a broad frequency range centered around delta’, Learning & Memory, Cold Spring Harbor Lab, 16, 69–81.

Park, P., Sanderson, T.M., Amici, M., Choi, S.-L., Bortolotto, Z.A., Zhuo, M., Kaang, B.-K. & Collingridge, G.L. (2016) ‘Calcium-Permeable AMPA Receptors Mediate the Induction of the Protein Kinase A-Dependent Component of Long-Term Potentiation in the Hippocampus’, Journal of Neuroscience, Society for Neuroscience, 36, 622–631.
