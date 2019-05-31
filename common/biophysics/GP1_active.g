//actpars - set conductance densities for the spiking dendrite model
//    Any param value can be overwritten subsequently.

//Voltage-gated ion channel densities
float G_Na_fast_GP      = 500 // changed from 350 12/10/2008      
float G_Na_slow_GP      = 1.015
float G_Kv3_GP          = 40    // changed from 11.25 12/13/2008
float G_Kv2_GP          = 20    // changed from 1 12/13/2008     
float G_Kv4_fast_GP     = 20 
float G_Kv4_slow_GP     = {G_Kv4_fast_GP}*1.5
float G_KCNQ_GP         = 2
float G_Ca_HVA_GP       = 0.3    
float G_K_ahp_GP        = 4
float G_h_HCN_GP        = 0.2
float G_h_HCN2_GP       = {G_h_HCN_GP}*2.5

//Multipliers for conductance densities
float G_mult             = 1
float G_mult_NaF_dend    = 1
float G_mult_NaP_dend    = 1.5
float G_mult_Kdr_dend      = 1
float G_mult_KA_dend       = 2
float G_mult_KCNQ_dend    = 1
float G_mult_SK_dend    = 0.03
float G_mult_Ca_dend     = 1
float G_mult_HCN_dend    = 1
float G_mult_NaF_soma    = 1
float G_mult_NaP_soma    = 1
float G_mult_Kdr_soma   = 1 
float G_mult_KA_soma    = 1
float G_mult_KCNQ_soma    = 1
float G_mult_SK_soma    = 16
float G_mult_Ca_soma    = 1
float G_mult_NaF_axon    = 3.5
float G_mult_NaP_axon    = 3.5
float G_mult_Kdr_axon     = 5    
float G_mult_KA_axon     = 5    
float G_mult_KCNQ_axon    = 5    

