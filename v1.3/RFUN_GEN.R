#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.3
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories
#' Prior knowledge from modern samples (gene polymorphism)

#' Genotype frequency data

#' R functions

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
sourceCpp("./Code/Code v1.0/Code 1L/Code v1.3/CFUN_GEN.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation of the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation of the simulated mutant allele frequency trajectory ended

#' Standard version
simulateWFM <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen) {
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
    WFM <- simulateWFM_arma(fts_mat, pop_siz, int_frq, int_gen, evt_gen)
    mut_frq_pth_pre_evt <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth_pre_evt <- as.matrix(WFM$gen_frq_pth)

    fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz, mut_frq_pth_pre_evt[, ncol(mut_frq_pth_pre_evt)], evt_gen, lst_gen)
    mut_frq_pth_pst_evt <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth_pst_evt <- as.matrix(WFM$gen_frq_pth)

    mut_frq_pth <- cbind(mut_frq_pth_pre_evt, mut_frq_pth_pst_evt[, -1])
    gen_frq_pth <- cbind(gen_frq_pth_pre_evt, gen_frq_pth_pst_evt[, -1])
  }

  return(list(mut_frq_pth = mut_frq_pth,
              gen_frq_pth = gen_frq_pth))
}
#' Compiled version
cmpsimulateWFM <- cmpfun(simulateWFM)

########################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param int_gen the generation of the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation of the simulated mutant allele frequency trajectory ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateWFD <- function(sel_cof, dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
  if (evt_gen >= lst_gen) {
    mut_frq_pth <- simulateWFD_arma(sel_cof[1], dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
    mut_frq_pth <- as.vector(mut_frq_pth)
  } else if (evt_gen < int_gen) {
    mut_frq_pth <- simulateWFD_arma(sel_cof[2], dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num)
    mut_frq_pth <- as.vector(mut_frq_pth)
  } else {
    mut_frq_pth_pre_evt <- simulateWFD_arma(sel_cof[1], dom_par, pop_siz, ref_siz, int_frq, int_gen, evt_gen, ptn_num)
    mut_frq_pth_pre_evt <- as.vector(mut_frq_pth_pre_evt)

    mut_frq_pth_pst_evt <- simulateWFD_arma(sel_cof[2], dom_par, pop_siz, ref_siz, tail(mut_frq_pth_pre_evt, n = 1), evt_gen, lst_gen, ptn_num)
    mut_frq_pth_pst_evt <- as.vector(mut_frq_pth_pst_evt)

    mut_frq_pth <- cbind(mut_frq_pth_pre_evt, mut_frq_pth_pst_evt[, -1])
  }

  if (dat_aug == FALSE) {
    return(mut_frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(mut_frq_pth)
  }
}
#' Compiled version
cmpsimulateWFD <- cmpfun(simulateWFD)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial mutant allele frequency of the population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # generate the population mutant allele frequency trajectory
  if (model == "WFM") {
    WFM <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_con, evt_gen, int_gen, lst_gen)
    pop_mut_frq <- as.vector(WFM$mut_frq_pth)
    pop_gen_frq <- as.matrix(WFM$gen_frq_pth)
  }
  if (model == "WFD") {
    pop_mut_frq <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_con, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
    pop_mut_frq <- as.vector(pop_mut_frq)
    pop_gen_frq <- matrix(NA, nrow = 10, ncol = ncol(pop_mut_frq))
    fts_mat <- calculateFitnessMat_arma(sel_cof, dom_par)
    for (k in 1:ncol(pop_mut_frq)) {
      ale_frq <- c(1 - pop_mut_frq[k], pop_mut_frq[k])
      gen_frq <- fts_mat * (ale_frq %*% t(ale_frq)) / sum(fts_mat * (ale_frq %*% t(ale_frq)))
      gen_frq[lower.tri(gen_frq, diag = FALSE)] <- NA
      pop_gen_frq[, k] <- discard(as.vector(2 * gen_frq - diag(diag(gen_frq), nrow = 4, ncol = 4)), is.na)
    }
    pop_gen_frq <- as.matrix(pop_gen_frq)
  }

  # generate the sample allele counts at all sampling time points
  smp_gen_cnt <- matrix(NA, nrow = 3, ncol = length(smp_gen))
  smp_gen_frq <- matrix(NA, nrow = 3, ncol = length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_gen_cnt[, k] <- rmultinom(1, size = smp_siz[k], prob = pop_gen_frq[, smp_gen[k] - int_gen + 1])
    smp_gen_frq[, k] <- smp_gen_cnt[, k] / smp_siz[k]
  }

  return(list(smp_gen = smp_gen,
              smp_siz = smp_siz,
              smp_gen_cnt = smp_gen_cnt,
              smp_gen_frq = smp_gen_frq,
              pop_gen_frq = pop_gen_frq,
              pop_mut_frq = pop_mut_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################

#' Run the bootstrap particle filter (BPF) with the single-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  # combine the sampling time points and the event time point
  if (evt_gen %in% smp_gen) {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]

    if (which(smp_gen == evt_gen)[1] == which(smp_siz == 0)) {
      smp_gen[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_gen[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_siz[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_siz[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_cnt[, c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_cnt[, c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
    }
  } else {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]
  }
  smp_gen <- smp_gen - min(smp_gen)

  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)

  return(list(lik = BPF$lik,
              wght = BPF$wght,
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
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  # combine the sampling time points and the event time point
  if (evt_gen %in% smp_gen) {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]

    if (which(smp_gen == evt_gen)[1] == which(smp_siz == 0)) {
      smp_gen[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_gen[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_siz[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_siz[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_cnt[, c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_cnt[, c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
    }
  } else {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]
  }
  smp_gen <- smp_gen - min(smp_gen)

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)

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
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  # combine the sampling time points and the event time point
  if (evt_gen %in% smp_gen) {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]

    if (which(smp_gen == evt_gen)[1] == which(smp_siz == 0)) {
      smp_gen[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_gen[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_siz[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_siz[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_cnt[, c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_cnt[, c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
    }
  } else {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]
  }
  smp_gen <- smp_gen - min(smp_gen)

  # run the PMMH
  sel_cof_chn <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  sel_cof_chn <- as.matrix(sel_cof_chn)

  return(sel_cof_chn)
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

########################################

#' Run the adaptive particle marginal Metropolis-Hastings (AdaptPMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runAdaptPMMH <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto) {
  # combine the sampling time points and the event time point
  if (evt_gen %in% smp_gen) {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]

    if (which(smp_gen == evt_gen)[1] == which(smp_siz == 0)) {
      smp_gen[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_gen[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_siz[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_siz[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_cnt[, c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_cnt[, c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
    }
  } else {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]
  }
  smp_gen <- smp_gen - min(smp_gen)

  # run the PMMH
  sel_cof_chn <- runAdaptPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  sel_cof_chn <- as.matrix(sel_cof_chn)

  return(sel_cof_chn)
}
#' Compiled version
cmprunAdaptPMMH <- cmpfun(runAdaptPMMH)

########################################

#' Run the Bayesian Procedure for the inference of the selection coefficient
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the genotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, ...) {
  # combine the sampling time points and the event time point
  if (evt_gen %in% smp_gen) {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]

    if (which(smp_gen == evt_gen)[1] == which(smp_siz == 0)) {
      smp_gen[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_gen[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_siz[c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_siz[c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
      smp_cnt[, c(which(smp_gen == evt_gen)[1], which(smp_gen == evt_gen)[2])] <- smp_cnt[, c(which(smp_gen == evt_gen)[2], which(smp_gen == evt_gen)[1])]
    }
  } else {
    smp_gen <- append(smp_gen, evt_gen)
    smp_siz <- append(smp_siz, 0)
    smp_cnt <- cbind(smp_cnt, matrix(0, nrow = nrow(smp_cnt), ncol = 1))

    odr_ind <- sort(smp_gen, decreasing = FALSE, index.return = TRUE)$ix
    smp_gen <- smp_gen[odr_ind]
    smp_siz <- smp_siz[odr_ind]
    smp_cnt <- smp_cnt[, odr_ind]
  }
  smp_gen <- smp_gen - min(smp_gen)

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    sel_cof_chn <- runAdaptPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    sel_cof_chn <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  }
  sel_cof_chn <- as.matrix(sel_cof_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
  sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]

  # MMSE estimate for selection coefficients
  sel_cof_est <- rowMeans(sel_cof_chn)

  # 95% HPD interval for selection coefficients
  sel_cof_hpd <- matrix(NA, nrow = 2, ncol = 2)
  sel_cof_hpd[1, ] <- HPDinterval(as.mcmc(sel_cof_chn[1, ]), prob = 0.95)
  sel_cof_hpd[2, ] <- HPDinterval(as.mcmc(sel_cof_chn[2, ]), prob = 0.95)

  # calculate the change in the selection coefficient before and after the event
  dif_sel_chn <- sel_cof_chn[2, ] - sel_cof_chn[1, ]

  # MMSE estimate for the change in the selection coefficient before and after the event
  dif_sel_est <- mean(dif_sel_chn)

  # 95% HPD intervals for the change in the selection coefficient before and after the event
  dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              dif_sel_est = dif_sel_est,
              dif_sel_hpd = dif_sel_hpd,
              dif_sel_chn = dif_sel_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
