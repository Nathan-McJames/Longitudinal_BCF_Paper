# Longitudinal BCF Paper

Materials for replicating the results in "Bayesian Causal Forests for Longitudinal Data: Assessing the Impact of Part-Time Work on Growth in High School Mathematics Achievement".

The preprint of this paper can be found at: https://arxiv.org/pdf/2407.11927.

The data used in this study is available for download at https://nces.ed.gov/datalab/onlinecodebook.

Further information on the dataset (and others organised by the NCES) can be found at https://nces.ed.gov/surveys/hsls09/. 

![alt text](https://github.com/Nathan-McJames/Longitudinal_BCF_Paper/blob/main/Pictures/dgp1_figure.svg?raw=true)

![alt text](https://github.com/Nathan-McJames/Longitudinal_BCF_Paper/blob/main/Pictures/dgp2_figure.svg?raw=true)

![alt text](https://github.com/Nathan-McJames/Longitudinal_BCF_Paper/blob/main/Pictures/growth_plot.svg?raw=true)

![alt text](https://github.com/Nathan-McJames/Longitudinal_BCF_Paper/blob/main/Pictures/treat_plot.svg?raw=true)

### Code Scripts and Description:

BCF_SIM_Github.R - Code for running the DGP1 simulation study.

LBCF_2_Test_Github.cpp - RCPP code that implements the two wave version of LBCF for use in DGP1.

<br/>

GEST_SIM_Github.R - Code for running the DGP1 simulation study.

LBCF_3_Test_Github.cpp - RCPP code that implements the three wave version of LBCF for use in DGP2.

<br/>

LBCF_Paper_1_Github.R - Code for running the main analysis of the paper, fitting LBCF to the HSLS data.

LBCF_HSLS_Github.cpp - RCPP code that implements the missing data version of LBCF used in the main analysis of the HSLS data.
