#' @title Estimating selection coefficients and testing their changes from ancient DNA data
#' @author Xiaoyang Dai, Mark Beaumont, Feng Yu, Zhangyi He

#' version 1.1
#' Phenotypes controlled by a single gene
#' Constant natural selection and non-constant demographic histories

#' Allele frequency data

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2019/HE2021-WFM-2L-DiffusApprox-PMMH1-MolEcol")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

#install.packages("emdbook")
library("emdbook")

# call R functions
source("./Code/Code v1.0/Code 1L/Code v1.1/RFUN_ALE.R")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the generation of the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation of the simulated mutant allele frequency trajectory ended

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFM: the mutant allele frequency trajectory")

########################################

#' Simulate the mutant allele frequency trajectory according to the single-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the generation of the simulated mutant allele frequency trajectory started
#' @param lst_gen the generation of the simulated mutant allele frequency trajectory ended
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param dat_aug = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth, type = 'l', lwd = 1.5,
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFD: the mutant allele frequency trajectory")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
ref_siz <- 1e+04
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

smp_WFM <- numeric(sim_num)
smp_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  smp_WFM[i] <- cmpsimulateWFM(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)[(lst_gen - int_gen) + 1]
  smp_WFD[i] <- cmpsimulateWFD(sel_cof, dom_par, pop_siz, ref_siz, int_frq, int_gen, lst_gen, ptn_num, dat_aug = FALSE)[(lst_gen - int_gen) + 1]
}

hist(smp_WFM, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
     xlim = c(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD)),
     xlab = "Allele frequency", main = paste("Histograms of the mutant allele frequency at generation", lst_gen))
hist(smp_WFD, breaks = seq(min(smp_WFM, smp_WFD), max(smp_WFM, smp_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param int_con the initial mutant allele frequency of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param ref_siz the reference size of the horse population
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model
model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
smp_frq <- sim_HMM_WFM$smp_frq
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq, type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFM-HMM: the mutant allele")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ref_siz <- 1e+04
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, smp_gen, smp_siz, ref_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_cnt <- sim_HMM_WFD$smp_cnt
smp_frq <- sim_HMM_WFD$smp_frq
pop_frq <- sim_HMM_WFD$pop_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq, type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "WFD-HMM: the mutant allele")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 1
set.seed(test_seed)

model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- c(rep(1e+04, length.out = 201), rep(5e+03, length.out = 200), rep(1e+04, length.out = 100))
int_con <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_con, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
smp_frq <- sim_HMM_WFM$smp_frq
pop_frq <- sim_HMM_WFM$pop_frq

save(model, sel_cof, dom_par, pop_siz, int_con, smp_gen, smp_siz, smp_cnt, smp_frq, pop_frq,
     file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_frq, type = 'l', lwd = 1.5,
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)),
     xlab = "Generation", ylab = "Allele frequency",
     main = "Mutant allele")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the single-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
ref_siz <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+05

system.time(BPF <- cmprunBPF(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num))

save(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, BPF,
     file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BPF.rda")

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BPF_Likelihood.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l',
     xlab = "Number of particles", ylab = "Log likelihood",
     main = "Log likelihood through the bootstrap particle filter")
dev.off()

pop_frq_pre_resmp <- BPF$ale_frq_pre_resmp
pop_frq_pst_resmp <- BPF$ale_frq_pst_resmp

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BPF_Particle.pdf", width = 32, height = 18)
par(mfrow = c(3, 4), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_frq_pst_resmp[, k], breaks = seq(min(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), max(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_frq_pre_resmp[, k], breaks = seq(min(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), max(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), length.out = 50), plot = FALSE)
  hist(pop_frq_pst_resmp[, k], breaks = seq(min(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), max(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5),
       xlim = c(min(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k], smp_frq[k]), max(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k], smp_frq[k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)),
       xlab = "Allele frequency",
       main = paste("Mutant allele at generation", smp_gen[k]))
  hist(pop_frq_pre_resmp[, k], breaks = seq(min(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), max(pop_frq_pst_resmp[, k], pop_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[k], col = 'red', lty = 2, lwd = 2)
}
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
ref_siz <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_OptNum.rda")

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_OptNum.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2,
     xlab = "Particle number", ylab = "Log-likelihood standard deviation",
     main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e+00
dom_par
pop_siz
ref_siz <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04

system.time(sel_cof_chn <- cmprunPMMH(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

save(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, sel_cof_chn,
     file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_PMMH.rda")

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_PMMH.rda")

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_PMMH_Traceplot.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l',
     xlab = "Iteration", ylab = "Selection coefficient",
     main = "Trace plot for selection coefficient")
abline(h = sel_cof[1], col = 'red', lty = 2, lwd = 2)
dev.off()

# brn_num <- 1e+04
brn_num <- 0.5 * length(sel_cof_chn) + 1
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 5e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

sel_cof_est <- mean(sel_cof_chn)

sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)

pdf(file = "./Output/Output v2.0/Test v1.0/TEST_1L_ALE_PMMH_Posterior.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian Procedure for the inference of the selection coefficient
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the size of the horse population (non-constant)
#' @param ref_siz the reference size of the horse population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the horses drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the PMMH
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

set.seed(test_seed)

sel_cof <- 0e+00
dom_par
pop_siz
ref_siz <- 1e+04
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 2e+04
brn_num <- 1e+04
thn_num <- 5e+00

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num))

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_SimData.rda")

save(sel_cof, dom_par, pop_siz, ref_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, BayesianProcedure,
     file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BayesProc.rda")

load("./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BayesProc.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn

sel_cof_est <- BayesianProcedure$sel_cof_est

sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./Output/Output v1.0/TEST v1.1/TEST_1L_ALE_BayesProc_Posterior.pdf", width = 8, height = 6)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient",
     main = "Posterior for selection coefficient")
lines(density(sel_cof_chn), lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_est, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
