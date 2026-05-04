# Optimal-Design-AC-SLSE
This project provides the R code required to reproduce the numerical results presented in the paper "Optimal Designs for Active-Controlled Dose-Response Models with Asymmetric Errors."

# figure
Produces Figure 1 in the paper.

# AC-optimal design and AC-efficiency 

Tab1_ACoptACeff_M-M.R implements the calculation of AC-optimal designs and AC-efficiencies for the Michaelis-Menten (M-M) model under different values of t_1.

Table 3 contains the calculation of the AC-optimal design based on the SLSE for the four candidate models, as well as the AC-efficiency of the AC-optimal design based on the OLSE.

Tab5_ACeffmis_t1.R implements the calculation of AC-efficiency w.r.t. mis-specified t_1 in the M-M model.

Table 6 contains the calculation of AC-efficiency w.r.t. mis-specified parameters in the M-M model.

# simulation
Simulation_M-M.R produces the results of the M-M model in Section 5.1

# Dependencies
This project requires R (version ≥ 4.4.0) and the ICAOD package. https://cran.r-project.org/web/packages/ICAOD/index.html
