# 4th-Year-MMath-Project

This repository contains the main codes used for simulations and the equation-free method described in my MMath project 'Investigating the potential of early warning signals in disease elimination'. 

The file 'gillespie_SIS' contains the Gillespie simulation of the SIS model.
Files beginning 'project_gillespie_...' followed by a data type generate the main results (approximated drift function and potential surface) for that data type. 
The file 'Error_estimation_prevalence' calculates the error in the approximated drift function and potential surface for prevalence against the analytic results as a function of the number of realisations of the simulation. 
The file 'poly_estimate' calculates the coefficients a polynomial (cubic) estimate of the bottom 10% of the potential surface. 
