README: Preregistered Code for Modelling Reaction Time and Confidence Distributions in Decision Making
====
Find the [preregistration on OSF](https://mfr.de-1.osf.io/render?url=https://osf.io/x548k/?direct%26mode=render%26action=download%26mode=render)


## Structure:

* dynWEV-source package file (.tar.gz)
* functions folder for small helper functions
* script R-file for the actual analyses, including:
  * Reading and Aggregating Data
  * Fitting model parameters and predict confidence and RT distributions
  * Compare model fits with BIC and AIC
  * Produce plots to visualise model fits

## Usage:

* Start R with package file in working directory `install.packages("dynWEV_0.0.tar.gz", repos = NULL, type = "source")`

* Source Script_Comparison_SeqSamplingConfidenceModels.R or run in RStudio
