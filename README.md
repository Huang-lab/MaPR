# **Master regulators governing protein abundance\**


**Functions**

This package has two functions: **Network** and **MaPR**.
The **Network** function constructs a regulatory network of protein translation independent of transcriptional regulation by simultaneously using transcriptomic data and proteomic data of the same cohort. <br /><br />
The **MaPR** function predicts master protein abundance regulators based on the permutation of the constructed protein translation regulatory network generated by the **Network** function.



**Installation**

The package can be installed directly from GitHub by typing the following in an R console:

if(!require("devtools")) install.packages("devtools")
devtools::install_github("https://github.com/Huang-lab/MaPR")
library(MaPR)
