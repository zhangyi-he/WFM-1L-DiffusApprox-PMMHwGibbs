#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Wenyang Lyu, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.6
#' Phenotypes controlled by a single gene
#' Non-constant natural selection and non-constant demographic histories

#' Integrate prior knowledge from modern samples (gene polymorphism)

#' Input: genotype likelihoods
#' Output: posteriors for the selection coefficient, the genotype frequency trajectories of the population and the genotypes of the sample

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
sourceCpp("./Code/Code v1.0/Code v1.6/CFUN.cpp")

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
    WFM <- simulateWFM_arma(fts_mat, pop_siz, int_frq, int_gen, evt_gen)
    mut_frq_pth_pre_evt <- as.vector(WFM$mut_frq_pth)
    gen_frq_pth_pre_evt <- as.matrix(WFM$gen_frq_pth)

    fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
    WFM <- simulateWFM_arma(fts_mat, pop_siz, mut_frq_pth_pre_evt[length(mut_frq_pth_pre_evt)], evt_gen, lst_gen)
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
simulateWFD <- function(sel_cof, dom_par, pop_siz, ref_siz, int_frq, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = TRUE) {
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

    mut_frq_pth <- append(mut_frq_pth_pre_evt, mut_frq_pth_pst_evt[-1])
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
#' @param mis_rat the rate of the missing allele observed in the sample
#' @param thr_val the threshold for genotype calling
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, pop_siz, int_con, evt_gen, smp_gen, smp_siz, mis_rat, thr_val, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)

  # generate the population genotype frequency trajectories
  los_cnt <- 0
  while (TRUE) {
    if (los_cnt > 1e+03) {
      message("Warning message: mutant allele lost completely from the population in all 1000 replicates.")
      break
    } else {
      los_cnt <- los_cnt + 1
    }

    if (model == "WFM") {
      WFM <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_con, evt_gen, int_gen, lst_gen)
      mut_frq <- as.vector(WFM$mut_frq_pth)
      gen_frq <- as.matrix(WFM$gen_frq_pth)
    }
    if (model == "WFD") {
      mut_frq <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_con, evt_gen, int_gen, lst_gen, ptn_num, dat_aug = FALSE)
      mut_frq <- as.vector(mut_frq)
      gen_frq <- matrix(NA, nrow = 3, ncol = length(mut_frq))
      if (evt_gen >= lst_gen) {
        fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par)
        for (k in 1:length(mut_frq)) {
          pop_ale_frq <- c(1 - mut_frq[k], mut_frq[k])
          pop_gen_frq <- fts_mat * (pop_ale_frq %*% t(pop_ale_frq)) / sum(fts_mat * (pop_ale_frq %*% t(pop_ale_frq)))
          pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
          gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 2, ncol = 2)), is.na)
        }
      } else if (evt_gen < int_gen) {
        fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
        for (k in 1:length(mut_frq)) {
          pop_ale_frq <- c(1 - mut_frq[k], mut_frq[k])
          pop_gen_frq <- fts_mat * (pop_ale_frq %*% t(pop_ale_frq)) / sum(fts_mat * (pop_ale_frq %*% t(pop_ale_frq)))
          pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
          gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 2, ncol = 2)), is.na)
        }
      } else {
        fts_mat <- calculateFitnessMat_arma(sel_cof[1], dom_par)
        for (k in 1:(evt_gen - int_gen + 1)) {
          pop_ale_frq <- c(1 - mut_frq[k], mut_frq[k])
          pop_gen_frq <- fts_mat * (pop_ale_frq %*% t(pop_ale_frq)) / sum(fts_mat * (pop_ale_frq %*% t(pop_ale_frq)))
          pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
          gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 2, ncol = 2)), is.na)
        }

        fts_mat <- calculateFitnessMat_arma(sel_cof[2], dom_par)
        for (k in (evt_gen - int_gen + 2):length(mut_frq)) {
          pop_ale_frq <- c(1 - mut_frq[k], mut_frq[k])
          pop_gen_frq <- fts_mat * (pop_ale_frq %*% t(pop_ale_frq)) / sum(fts_mat * (pop_ale_frq %*% t(pop_ale_frq)))
          pop_gen_frq[lower.tri(pop_gen_frq, diag = FALSE)] <- NA
          gen_frq[, k] <- discard(as.vector(2 * pop_gen_frq - diag(diag(pop_gen_frq), nrow = 2, ncol = 2)), is.na)
        }
      }
      gen_frq <- as.matrix(gen_frq)
    }

    if (tail(mut_frq, 1) > 0 && tail(mut_frq, 1) < 1) {
      break
    }
  }

  # generate the sample genotype counts at all sampling time points
  raw_smp <- NULL
  for (k in 1:length(smp_gen)) {
    raw_smp <- cbind(raw_smp, rbind(rep(smp_gen[k], times = smp_siz[k]), rmultinom(smp_siz[k], size = 1, prob = gen_frq[, smp_gen[k] - int_gen + 1])))
  }

  imp_smp <- raw_smp
  imp_smp <- as.data.frame(t(imp_smp))
  rownames(imp_smp) <- NULL
  colnames(imp_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

  mis_smp <- rmultinom(ncol(raw_smp), size = 1, prob = c((1 - mis_rat)^2, 2 * (1 - mis_rat) * mis_rat, mis_rat^2))
  for (i in 1:ncol(mis_smp)) {
    if (mis_smp[1, i] == 1) {
      idx <- c(which(raw_smp[2:4, i] == 1), which(raw_smp[2:4, i] != 1)) + 1
      raw_smp[idx[1], i] <- runif(n = 1, min = thr_val, max = 1)
      raw_smp[idx[2], i] <- runif(n = 1, min = 0, max = 1 - raw_smp[idx[1], i])
      raw_smp[idx[3], i] <- 1 - raw_smp[idx[1], i] - raw_smp[idx[2], i]
    }
    if (mis_smp[2, i] == 1) {
      if (raw_smp[2, i] == 1) {
        raw_smp[2, i] <- runif(n = 1, min = 0.5, max = thr_val)
        raw_smp[3, i] <- 1 - raw_smp[2, i]
      } else if (raw_smp[4, i] == 1) {
        raw_smp[4, i] <- runif(n = 1, min = 0.5, max = thr_val)
        raw_smp[3, i] <- 1 - raw_smp[2, i]
      } else {
        raw_smp[3, i] <- runif(n = 1, min = 0.5, max = thr_val)
        raw_smp[2, i] <- runif(n = 1, min = 0, max = 1 - raw_smp[3, i])
        raw_smp[4, i] <- 1 - raw_smp[2, i] - raw_smp[3, i]
      }
    }
    if (mis_smp[3, i] == 1) {
      idx <- c(which(raw_smp[2:4, i] == 1), which(raw_smp[2:4, i] != 1)) + 1
      raw_smp[idx[1], i] <- runif(n = 1, min = 0.5, max = thr_val)
      raw_smp[idx[2], i] <- runif(n = 1, min = 0, max = 1 - raw_smp[idx[1], i])
      raw_smp[idx[3], i] <- 1 - raw_smp[idx[1], i] - raw_smp[idx[2], i]
    }
  }
  raw_smp <- as.data.frame(t(raw_smp))
  rownames(raw_smp) <- NULL
  colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

  return(list(raw_smp = raw_smp,
              imp_smp = imp_smp,
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
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param frq_pth the mutant allele frequency of the population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 4) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num)

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
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param evt_gen the generation that the event of interest occurred
#' @param frq_pth the mutant allele frequency of the population
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 4) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # calculate the optimal particle number
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, frq_pth, raw_smp, ptn_num, pcl_num, gap_num)

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
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

#' Standard version
runPMMH <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 4) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)

  return(list(sel_cof_chn = as.matrix(PMMH$sel_cof_chn),
              frq_pth_chn = as.matrix(PMMH$frq_pth_chn),
              imp_smp_chn = as.array(PMMH$imp_smp_chn)))
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
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runAdaptPMMH <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 4) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  # run the adaptive PMMH
  PMMH <- runAdaptPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)

  return(list(sel_cof_chn = as.matrix(PMMH$sel_cof_chn),
              frq_pth_chn = as.matrix(PMMH$frq_pth_chn),
              imp_smp_chn = as.array(PMMH$imp_smp_chn)))
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
#' @param raw_smp the sample of ancient horses (sampling times and individual genotypes)
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param adp_set = TRUE/FALSE (return the result with the adaptive setting or not)
#' @param stp_siz the step size sequence in the adaptive setting (decaying to zero)
#' @param apt_rto the target mean acceptance probability of the adaptive setting

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, brn_num, thn_num, adp_set, ...) {
  # preprocess the raw sample
  if (ncol(raw_smp) == 4) {
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  } else {
    raw_smp <- raw_smp[, -(2:3)]
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- c("generation", "A0A0", "A0A1", "A1A1")

    raw_smp <- raw_smp[order(raw_smp$generation), ]
    evt_gen <- evt_gen - min(raw_smp$generation)
    raw_smp$generation <- raw_smp$generation - min(raw_smp$generation)
    rownames(raw_smp) <- NULL
    colnames(raw_smp) <- NULL
    raw_smp <- t(as.matrix(raw_smp))
  }

  if (adp_set == TRUE) {
    # run the adaptive PMMH
    PMMH <- runAdaptPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num, stp_siz, apt_rto)
  } else {
    # run the PMMH
    PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, ref_siz, evt_gen, raw_smp, ptn_num, pcl_num, itn_num)
  }
  sel_cof_chn <- as.matrix(PMMH$sel_cof_chn)
  frq_pth_chn <- as.matrix(PMMH$frq_pth_chn)
  imp_smp_chn <- as.array(PMMH$imp_smp_chn)

  # burn-in and thinning
  sel_cof_chn <- sel_cof_chn[, brn_num:dim(sel_cof_chn)[2]]
  sel_cof_chn <- sel_cof_chn[, (1:round(dim(sel_cof_chn)[2] / thn_num)) * thn_num]
  frq_pth_chn <- frq_pth_chn[brn_num:dim(frq_pth_chn)[1], ]
  frq_pth_chn <- frq_pth_chn[(1:round(dim(frq_pth_chn)[1] / thn_num)) * thn_num, ]
  imp_smp_chn <- imp_smp_chn[, , brn_num:dim(imp_smp_chn)[3]]
  imp_smp_chn <- imp_smp_chn[, , (1:round(dim(imp_smp_chn)[3] / thn_num)) * thn_num]

  # MMSE estimate for selection coefficients and mutant allele frequencies
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
  dif_sel_est <- mean(dif_sel_chn)

  # 95% HPD intervals for the change in the selection coefficient before and after the event
  dif_sel_hpd <- HPDinterval(as.mcmc(dif_sel_chn), prob = 0.95)

  # posterior probability for genotypes
  imp_smp_est <- matrix(NA, nrow = ncol(raw_smp), ncol = 3)
  for (i in 1:nrow(imp_smp_est)) {
    imp_smp_est[i, ] <- rowSums(imp_smp_chn[, i, ]) / dim(imp_smp_chn)[3]
  }
  imp_smp_est <- cbind(raw_smp[1, ], imp_smp_est)
  imp_smp_est <- as.data.frame(imp_smp_est)
  rownames(imp_smp_est) <- NULL
  colnames(imp_smp_est) <- c("generation", "A0A0", "A0A1", "A1A1")

  return(list(sel_cof_est = sel_cof_est,
              sel_cof_hpd = sel_cof_hpd,
              sel_cof_chn = sel_cof_chn,
              dif_sel_est = dif_sel_est,
              dif_sel_hpd = dif_sel_hpd,
              dif_sel_chn = dif_sel_chn,
              frq_pth_est = frq_pth_est,
              frq_pth_hpd = frq_pth_hpd,
              frq_pth_chn = frq_pth_chn,
              imp_smp_est = imp_smp_est,
              imp_smp_chn = imp_smp_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
