# Estimating temporally variable selection intensity from ancient DNA data II
The source code is implemented for an MCMC-based method for inferring temporally variable selection from ancient DNA data (in the format of genotype likelihoods), and the manuscript has been submitted to Molecular Biology and Evolution, available at https://doi.org/10.1101/2023.07.10.548348. Except for the PMMH-within-Gibbs procedure described in the manuscript, we also provide an adaptive version that improves the convergence of the chain.

[Code/FwdSimu](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-PMMHwGibbs/tree/main/Code/FwdSimu/) includes the source code implemented for the version with a full forward-in-time simulation of the Wright-Fisher diffusion.

[Code/BwdSimu](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-PMMHwGibbs/tree/main/Code/BwdSimu/) includes the source code implemented for the version with a full backward-in-time simulation of the Wright-Fisher diffusion.

[Code/MixSimu](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-PMMHwGibbs/tree/main/Code/MixSimu/) includes the source code implemented for the version with a mix of the forward- and backward-in-time simulations of the Wright-Fisher diffusion.

[Data](https://github.com/zhangyi-he/WFM-1L-DiffusApprox-PMMHwGibbs/tree/main/Data) includes the ancient horse samples genotyped at the loci for coat colouration.
