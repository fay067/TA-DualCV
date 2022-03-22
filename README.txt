This folder contains the implementation of the following paper:
TA-DualCS: Time-Aware Dual-Cross-Sectional Imputation for Multivariate Time Series Clinical Data (paper ID: 4628)

-----------------------------Files Directory-----------------------------
code_paper4628
|
|--code files
|--data                                 * Put the downloaded datasets here. We provide example of MIMIC-III dataset.
|    |--mimic
|         |
|         |--data_groundtruth
|         |
|         |--data_with_missing
| 
| 
|--result                             * The imputation results and prediction results are here.
|    |--mimic
|         |--mimic_24h_imputation_result
|         |
|         |--mimic_24h_prediction



-----------------------------Code-----------------------------
- For imputation, install the following dependencies:
  - hash
  - doParallel
  - foreach
  - abind
  - MICE
  - GPfit
- For prediction, install the following dependencies under Python3:
  - sys, math, numpy=1.18.1, pandas=1.0.1, collections
  - tensorflow=2.4.1
  - sklearn=0.22.1
 

1. Generate imputation results.
Rscript TADualCS.R

2. Septic shock 24-hour early prediction.
python septic_prediction.py [dataset] [missing indicator:1 for with MI; 0 for without MI]

An example:
python septic_prediction.py mimic 1
