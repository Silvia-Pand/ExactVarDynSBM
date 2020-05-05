This .zip file contains the R functions used for the paper
"An exact algorithm for time-dependent variational inference for the dynamic stochastic block model" by
- F.Bartolucci (University of Perugia, IT)
- S.Pandolfi (University of Perugia, IT)  

est_var_dyn_exact.R 		estimate the dynamic SBMs using the exact maximization 					algorithm for variational inference

est_var_dyn_approximate.R	estimate the dynamic SBMs using the variational 					approximation proposed by Matias and Miele (2017)
 	
example.R 			example file that loads the workspace file 						"example_data.RData" and fits the dynamic SBM using both 				est_var_dyn_exact.R and est_var_dyn_approximate.R  

example_data.RData  		workspace file containing a simulated dynamic network

draw_sn_dyn.R 			simulate dynamic SBMs

example_sim.R			example file that calls draw_sn_dyn.R and then fits the 				dynamic SBM using est_var_dyn_exact.R  
