% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fm_fit_model.R
\name{fm_fit_model}
\alias{fm_fit_model}
\title{Compile model and fit to data}
\usage{
fm_fit_model(
  fm_data_inputs,
  SE,
  data_source,
  data_obs,
  samp_process,
  b_constraint,
  cov_samp_process,
  biomass_temporal,
  sampling_temporal,
  lf_link,
  ref_data = c("com", "sci"),
  EM = c("est_b", "fix_b"),
  month_ref,
  compute_sd = FALSE,
  Version = "",
  seed = 29510
)
}
\arguments{
\item{fm_data_inputs}{list Named list obtained from `fm_load_data()`}

\item{SE}{Numeric. Apply bias correction and compute standard errors for key quantities of interest (spatial predictions and total biomass). N.b. can be very time consuming. - 0: no, 1: yes}

\item{data_source}{Numeric. Specify the data sources used to fit the model.    1 : integrated model (scientific + commercial data) ; 2 : scientific model ; 3 : commercial model}

\item{data_obs}{Numeric. Observation model for biomass.    1 : zero-inflated gamma ; 2 : zero-inflated lognormal ; 3 : lognormal}

\item{samp_process}{Numeric. Include sampling process or not.    0 : no sampling process (faster) ; 1 : inhomogeneous Poisson point process}

\item{b_constraint}{Numeric. Put constraint on b parameters.    1 : b are positives ; 2 : no constraints}

\item{cov_samp_process}{Numeric. Covariates in the sampling process.    0 : none ; 1 : covariates in the sampling process}

\item{biomass_temporal}{Numeric. Temporal correlation in biomass.    0 : no ; 1 : first-order autoregressive model (AR1)}

\item{sampling_temporal}{Numeric. Temporal correlation in sampling process.    0 : no ; 1 : first-order autoregressive model (AR1)}

\item{lf_link}{Numeric. Link function of the latent field.    0 : log (biomass data) ; 1 : logit (presence absence data)}

\item{ref_data}{character. Reference data for estimating the interecept of the latent field.    com : commercial (default) ; sci : scientific}

\item{EM}{character. Reference data for estimating the intercept of the latent field.    est_b : b is estimated ; fix_b : b is fixed}

\item{month_ref}{numeric.    Reference month    1 to 12}

\item{compute_sd}{logical    run sdreport().    TRUE ; FALSE}

\item{Version}{version for fm_build_data_params_map}

\item{seed}{integer The seed controlling for random effect. Default is 29510.}
}
\value{
list A named list of all necessary outputs for generating graphs (step 3)
}
\description{
Compile model and fit to data
}
\examples{
# run part 2
fm_data_inputs <- readr::read_rds(
  system.file("examples", "part1_output_small.rds",
              package = "FishMap")
  )

fm_model_results <- fm_fit_model(fm_data_inputs = fm_data_inputs,
                                 SE = 1,
                                 data_source = 1,
                                 data_obs = 2,
                                 samp_process = 0,
                                 b_constraint = 2,
                                 cov_samp_process = 0,
                                 biomass_temporal = 1,
                                 sampling_temporal = 0,
                                 lf_link = 0,
                                 ref_data = "com",
                                 EM = "est_b",
                                 month_ref = 1)
}
