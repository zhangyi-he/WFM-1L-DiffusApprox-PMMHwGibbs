#' @title Estimating selection coefficients and testing their changes from ancient DNA data II
#' @author Zhangyi He, Wenyang Lyu, Xiaoyang Dai, Mark Beaumont, Feng Yu

#' version 1.0 (backward-in-time simulation)
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient, the genotype frequency trajectories of the population and the genotypes of the sample

#' R functions

#install.packages("gtools")
library("gtools")

#install.packages("purrr")
library("purrr")

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./Code/Code v1.0/Code BkdSimu v1.0/CFUN.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended

#' Standard version
simulateWFM <- function(sel_cof, dom_par, pop_siz, int_frq, evt_gen, int_gen, lst_gen) {
  if (evt_gen >= lst_gen) {
    fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz, int_frq, int_gen, lst_gen)
    mut_frq_pth <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth <- as.matrix(WFM$gen_frq_pth)
  } else if (evt_gen < int_gen) {
    fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz, int_frq, int_gen, lst_gen)
    mut_frq_pth <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth <- as.matrix(WFM$gen_frq_pth)
  } else {
    fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz[1:(evt_gen - int_gen + 1)], int_frq, int_gen, evt_gen)
    mut_frq_pth_pre_evt <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth_pre_evt <- as.matrix(WFM$gen_frq_pth)

    fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz[(evt_gen - int_gen + 1):(lst_gen - int_gen + 1)], tail(mut_frq_pth_pre_evt, n = 1), evt_gen, lst_gen)
    mut_frq_pth_pst_evt <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth_pst_evt <- as.matrix(WFM$gen_frq_pth)

    mut_frq_pth <- append(mut_frq_pth_pre_evt, mut_frq_pth_pst_evt[-1])
    gen_frq_pth <- cbind(gen_frq_pth_pre_evt, gen_frq_pth_pst_evt[, -1])
  }

  return(list(mut_frq_pth = mut_frq_pth,
              gen_frq_pth = gen_frq_pth))
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the stochastic differential equation resulting from the single-locus Wright-Fisher diffusion with selection backward in time
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended
#' @param ptn_num
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, dom_par, pop_siz, ref_siz, lst_frq, evt_gen, int_gen, lst_gen, ptn_num, unif_grd, dat_aug = TRUE) {
  coord_x = calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)

  if (evt_gen >= lst_gen) {
    WFD <- simulateWFD_arma(sel_cof[1], dom_par, pop_siz, ref_siz, lst_frq, int_gen, lst_gen, ptn_num, coord_x)
    mut_frq_pth <- as.vector(WFD$mut_frq_pth)
    log_w <- as.double(WFD$log_w)
  } else if (evt_gen < int_gen) {
    WFD <- simulateWFD_arma(sel_cof[2], dom_par, pop_siz, ref_siz, lst_frq, int_gen, lst_gen, ptn_num, coord_x)
    mut_frq_pth <- as.vector(WFD$mut_frq_pth)
    log_w <- as.double(WFD$log_w)
  } else {
    WFD <- simulateWFD_arma(sel_cof[2], dom_par, pop_siz[(evt_gen - int_gen + 1):(lst_gen - int_gen + 1)], ref_siz, lst_frq, evt_gen, lst_gen, ptn_num, coord_x)
    mut_frq_pth_pst_evt <- as.vector(WFD$mut_frq_pth)
    wght_pst_evt <- as.vector(WFD$wght)

    WFD <- simulateWFD_arma(sel_cof[1], dom_par, pop_siz[1:(evt_gen - int_gen + 1)], ref_siz, head(mut_frq_pth_pst_evt, n = 1), int_gen, evt_gen, ptn_num, coord_x)
    mut_frq_pth_pre_evt <- as.vector(WFD$mut_frq_pth)
    wght_pre_evt <- as.vector(WFD$wght)

    mut_frq_pth <- append(mut_frq_pth_pre_evt, mut_frq_pth_pst_evt[-1])
    wght <- append(wght_pre_evt[-1], wght_pst_evt)
  }

  if (dat_aug == FALSE) {
    return(list(mut_frq_pth = mut_frq_pth[(0:(lst_gen - int_gen)) * ptn_num[1] + 1],
                wght = wght[(0:(lst_gen - int_gen)) * ptn_num[1] + 1]))
  } else {
    return(list(mut_frq_pth = mut_frq_pth,
                wght = wght))
  }
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Solve the partial differential equation resulting from the single-locus Wright-Fisher diffusion with selection backward in time
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation that the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation that the simulated mutant allele frequency trajectory ended
#' @param ptn_num
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
solveWFD <- function(sel_cof, dom_par, pop_siz, ref_siz, lst_frq_dist, evt_gen, int_gen, lst_gen, ptn_num, unif_grd, dat_aug = TRUE) {
  coord_x = calculateGrid_arma(pop_siz, ref_siz, ptn_num, unif_grd)

  if (evt_gen >= lst_gen) {
    sol <- solveWFD_arma(sel_cof[1], dom_par, pop_siz, ref_siz, lst_frq_dist, int_gen, lst_gen, ptn_num, coord_x)
    sol <- as.matrix(sol)
  } else if (evt_gen < int_gen) {
    sol <- solveWFD_arma(sel_cof[2], dom_par, pop_siz, ref_siz, lst_frq_dist, int_gen, lst_gen, ptn_num, coord_x)
    sol <- as.matrix(sol)
  } else {
    sol <- solveWFD_arma(sel_cof[2], dom_par, pop_siz[(evt_gen - int_gen + 1):(lst_gen - int_gen + 1)], ref_siz, lst_frq_dist, evt_gen, lst_gen, ptn_num, coord_x)
    sol_pst_evt <- as.matrix(sol)

    sol <- solveWFD_arma(sel_cof[1], dom_par, pop_siz[1:(evt_gen - int_gen + 1)], ref_siz, sol_pst_evt[, 1], int_gen, evt_gen, ptn_num, coord_x)
    sol_pre_evt <- as.matrix(sol)

    sol <- cbind(sol_pre_evt, sol_pst_evt[, -1])
  }

  if (dat_aug == FALSE) {
    return(sol[, (0:(lst_gen - int_gen)) * ptn_num[1] + 1])
  } else {
    return(sol)
  }
}
#' Compiled version
cmpsolveWFD <- cmpfun(solveWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_lab the identifier of the sample assigned
#' @param smp_gen the generation of the sample drawn
#' @param smp_qua the quality of the sample tested
#' @param thr_val the threshold for genotype calling

#' Standard version
simulateHMM <- function(sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_lab, smp_gen, smp_qua, thr_val) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # generate the population genotype frequency trajectories
  WFM <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_con, evt_gen, int_gen, lst_gen)
  mut_frq <- as.vector(WFM$mut_frq_pth)
  gen_frq <- as.matrix(WFM$gen_frq_pth)

  # generate the sample genotype counts at all sampling time points
  tru_smp <- NULL
  for (k in 1:length(smp_gen)) {
    tru_smp <- cbind(tru_smp, rbind(rep(smp_gen[k], times = 1), rmultinom(1, size = 1, prob = gen_frq[, smp_gen[k] - int_gen + 1])))
  }
  tru_smp <- as.data.frame(t(tru_smp))
  rownames(tru_smp) <- smp_lab
  colnames(tru_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  tru_smp <- tru_smp[order(tru_smp$generation), ]

  raw_smp <- tru_smp
  sel_idx <- which(tru_smp$'A0A0' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c(smp_qua[1], (1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'A0A1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, smp_qua[1], (1 - smp_qua[1]) / 2) * smp_qua[2])
  sel_idx <- which(tru_smp$'A1A1' == 1)
  raw_smp[sel_idx, 2:4] <- rdirichlet(length(sel_idx), alpha = c((1 - smp_qua[1]) / 2, (1 - smp_qua[1]) / 2, smp_qua[1]) * smp_qua[2])

  cal_smp <- raw_smp
  mis_idx <- 1:nrow(raw_smp)
  sel_idx <- which(raw_smp$'A0A0' / raw_smp$'A0A1' >= thr_val & raw_smp$'A0A0' / raw_smp$'A1A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(1, 0, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'A0A1' / raw_smp$'A0A0' >= thr_val & raw_smp$'A0A1' / raw_smp$'A1A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 1, 0), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  sel_idx <- which(raw_smp$'A1A1' / raw_smp$'A0A0' >= thr_val & raw_smp$'A1A1' / raw_smp$'A0A1' >= thr_val)
  if (length(sel_idx) > 1) {
    cal_smp[sel_idx, 2:4] <- matrix(rep(c(0, 0, 1), times = length(sel_idx)), nrow = length(sel_idx), ncol = 3, byrow = TRUE)
    mis_idx <- mis_idx[-which(mis_idx %in% sel_idx)]
  }
  if (length(mis_idx) > 1) {
    cal_smp[mis_idx, 2:4] <- matrix(rep(c(NA, NA, NA), times = length(mis_idx)), nrow = length(mis_idx), ncol = 3, byrow = TRUE)
  }

  fil_smp <- cal_smp[complete.cases(cal_smp), ]

  cal_rat <- nrow(fil_smp) / nrow(cal_smp)
  # cal_rat
  err_rat <- sum(rowSums(cal_smp[, -(1:2)] == tru_smp[, -(1:2)]) == 1, na.rm = TRUE) / nrow(fil_smp)
  # err_rat

  return(list(tru_smp = tru_smp,
              raw_smp = raw_smp,
              cal_smp = cal_smp,
              fil_smp = fil_smp,
              cal_rat = cal_rat,
              err_rat = err_rat,
              gen_frq = gen_frq,
              mut_frq = mut_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################

#' Run the bootstrap particle filter (BPF) with the single-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mut_rat the mutation rate
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param frq_pth the mutant allele frequency trajectory of the population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num)

  return(list(lik = BPF$lik,
              wght = BPF$wght,
              mut_frq_pth = BPF$mut_frq_pth,
              mut_frq_pre_resmp = BPF$mut_frq_pre_resmp,
              mut_frq_pst_resmp = BPF$mut_frq_pst_resmp,
              gen_frq_pre_resmp = BPF$gen_frq_pre_resmp,
              gen_frq_pst_resmp = BPF$gen_frq_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mut_rat the mutation rate
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param frq_pth the mutant allele frequency trajectory of the population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num)

  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num),
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mut_rat the mutation rate
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num)

  return(list(sel_cof_chn = as.matrix(PMMH$sel_cof_chn),
              frq_pth_chn = as.matrix(PMMH$frq_pth_chn),
              cal_smp_chn = as.array(PMMH$cal_smp_chn)))
}
#' Compiled version
cmprunBackPMMH <- cmpfun(runPMMH)

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mut_rat the mutation rate
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runAdaptPMMH <- function(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  # run the adaptive PMMH
  PMMH <- runAdaptPMMH_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto)

  return(list(sel_cof_chn = as.matrix(PMMH$sel_cof_chn),
              frq_pth_chn = as.matrix(PMMH$frq_pth_chn),
              cal_smp_chn = as.array(PMMH$cal_smp_chn)))
}
#' Compiled version
cmprunAdaptBackPMMH <- cmpfun(runAdaptPMMH)

########################################

#' Run the Bayesian Procedure for the inference of the selection coefficient
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param mut_rat the mutation rate
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, brn_num, thn_num, grd_num, adp_set, ...) {
  # preprocess the raw sample
  raw_smp <- as.data.frame(raw_smp)
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")
  raw_smp <- raw_smp[order(raw_smp$generation), ]
  smp_lab <- rownames(raw_smp)
  smp_gen <- raw_smp$generation
  evt_gen <- evt_gen - min(raw_smp$generation)
  raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)

  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- NULL
  raw_smp <- t(as.matrix(raw_smp))

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    PMMH <- runAdaptPMMH_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    PMMH <- runPMMH_arma(sel_cof, dom_par, mut_rat, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, unif_grd, pcl_num, itn_num)
  }
  sel_cof_chn <- as.matrix(PMMH$sel_cof_chn)
  frq_pth_chn <- as.matrix(PMMH$frq_pth_chn)
  cal_smp_chn <- as.array(PMMH$cal_smp_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
  sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
  frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
  frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
  cal_smp_chn <- cal_smp_chn[, , brn_num:dim(cal_smp_chn)[3]]
  cal_smp_chn <- cal_smp_chn[, , (1:round(dim(cal_smp_chn)[3] / thn_num)) * thn_num]

  # MMSE estimates for selection coefficients and mutant allele frequencies
  # sel_cof_pdf <- kde2d(sel_cof_chn[1, ], sel_cof_chn[2, ], n = grd_num)
  # sel_cof_est <- c(sel_cof_pdf$x[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]], sel_cof_pdf$y[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]])
  sel_cof_est <- rowMeans(sel_cof_chn)
  frq_pth_est <- colMeans(frq_pth_chn)

  # 95% HPD interval for selection coefficients and mutant allele frequencies
  sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
  sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
  sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)
  frq_pth_hpd <- matrix(NA, nrow = 2, ncol = dim(frq_pth_chn)[2])
  for (i in 1:dim(frq_pth_chn)[2]) {
    frq_pth_hpd[, i] <- HPDinterval(as.mcmc(frq_pth_chn[, i]), prob = 0.95)
  }

  # calculate the change in the selection coefficient before and after the event
  dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

  # MMSE estimate for the change in the selection coefficient before and after the event
  # dif_sel_pdf <- density(dif_sel_chn, n = grd_num)
  # dif_sel_est <- dif_sel_pdf$x[which(dif_sel_pdf$y == max(dif_sel_pdf$y))]
  dif_sel_est <- mean(dif_sel_chn)

  # 95% HPD interval for the change in the selection coefficient before and after the event
  dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

  # posterior probability for genotypes
  cal_smp_est <- matrix(NA, nrow = ncol(raw_smp), ncol = 3)
  for (i in 1:nrow(cal_smp_est)) {
    cal_smp_est[i, ] <- rowSums(cal_smp_chn[, i, ]) / dim(cal_smp_chn)[3]
  }
  cal_smp_est <- as.data.frame(cal_smp_est)
  cal_smp_est <- cbind(smp_gen, cal_smp_est)
  rownames(cal_smp_est) <- smp_lab
  colnames(cal_smp_est) <- c("generation", "A0A0", "A0A1", "A1A1")

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              dif_sel_est = dif_sel_est,
              dif_sel_hpd = dif_sel_hpd,
              dif_sel_chn = dif_sel_chn,
              frq_pth_est = frq_pth_est,
              frq_pth_hpd = frq_pth_hpd,
              frq_pth_chn = frq_pth_chn,
              cal_smp_est = cal_smp_est,
              cal_smp_chn = cal_smp_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
