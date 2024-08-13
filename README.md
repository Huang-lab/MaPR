# **Master regulators governing protein abundance\**


### Introduction
This package has two functions: **Network** and **MaPR**.<br /><br />
The **Network** function constructs a regulatory network of protein translation independent of transcriptional regulation by simultaneously using transcriptomic data and proteomic data of the same cohort. <br /><br />
The **MaPR** function predicts master protein abundance regulators based on the permutation of the constructed protein translation regulatory network generated by the **Network** function.

<br /><br />
### Installation
The package can be installed directly from GitHub by typing the following in an `R` console:<br /><br />
```R
if(!require("devtools")) install.packages("devtools")

devtools::install_github("https://github.com/Huang-lab/MaPR")
library(MaPR)
```

<br /><br />
### Citation
If you used or adapted MaPR in your study, please cite our paper [].

<br /><br />
### Contact
If you run into issues, you can contact authors: zishan.wang{at}mssm.edu or kuan-lin.huang{at}mssm.edu.
