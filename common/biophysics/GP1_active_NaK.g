//Set conductance densities for the model with uniformly high dendritic gNaF
//	but only 2 ion channels: NaF and Kv2
//	Any param value can be overwritten subsequently.

//Voltage-gated ion channel densities
float G_Na_fast_GP  = 500	 
float G_Kv2_GP  	= 60	

//Multipliers for conductance densities
float G_mult 			= 1
float G_mult_NaF_dend   = 1
float G_mult_Kdr_dend  	= 1
float G_mult_NaF_soma   = 1
float G_mult_Kdr_soma   = 1
float G_mult_NaF_axon 	= 40
float G_mult_Kdr_axon 	= 40

float ELEAK_sd = -0.057
float ELEAK_ax = {ELEAK_sd}
float EREST_ACT = {ELEAK_sd}
