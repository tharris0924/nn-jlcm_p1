# Deep learning embedded latent class joint models of time-to-event and longitudinal data
The code in this repo showcases the method proposed by Harris, Nakhaeirad and Skhosana (2026). The main function is the `nnem.R` function. This function is used to compute a joint latent class model for time-to-event and longitudinal data using an adapted EM algorithm. This EM algorithm, uses a neural network to estimate the prior mixing proportions conditonal on a subjects covariates. This approach allows users to flexibly model the mixing proportions without relying on the assumptions of a logistic regression model, which is used in the traditional approach (see Proust-Lima et. al (2017)). 

Here is a display of the mixing proportions estimated as the algorithm runs for 20 iterations on simulation III for K=4 components:

![](https://github.com/tharris0924/nn-jlcm_p1/blob/main/Codes/mixing_proportions.gif)

This repo provides the code to be able to generate the data for the simulated data above in `gen_data_k4.R` as well as the data used in the application in the `gen_data_paquid.R` function.

## References

Harris, T. J. E., Nakhaei Rad, N., & Skhosana, S. B. (2026). Deep learning embedded latent class joint modelling of time-to-event and longitudinal data. Statistical Methods in Medical Research, 1–39. https://doi.org/10.1177/09622802261465334 <br>

Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of extended mixed models using latent classes and latent processes: The R package lcmm. Journal of Statistical Software, 78(2), 1–56. https://doi.org/10.18637/jss.v078.i02
