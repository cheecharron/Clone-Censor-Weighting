# Clone-Censor-Weighting
Application of clone censor weighting to longitudinal data for timing-of-exposures analyses

The 'simulate data.R' file can be used to simulate longitudinal data that contains data on 10,000
pseudo-patients, with 24 months of monthly data on time-varying weight-loss medication (WLM), 
weight, BMI, diabetes status, and change from initial weight. The dataset also contains non-time-
varying info on patients' height and group membership ("early", "late", or "never") based on 
timing of WLM.

The 'clone censor weighting - Cox.R' file is used to generate clone censor weights, using Cox
survival modeling. This code works by creating clones for each patient, one for each clone. 
Censoring is used when the clone violates group membership. Cox modeling is applied to each group
to estimate the probability of censoring at each time point. A numerator model (with time-varying
and time-invariant predictors) and denominator model (with time-invariant predictors only) are
specified for each group. The inverse probabilities of an event occurring by each time point
from the numerator model are divided by those from the denominator model. The datasets for each
group are then merged into one long dataset.

