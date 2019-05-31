//Passive membrane properties for GP1
float RA        = 1.74      // uniform
float CM        = 0.024     // all unmyelinated regions
float CM_my     = 0.00024   // myelinated axon segments.
float RM_sd     = 1.47      // soma
float RM_ax     = {RM_sd}   // unmyelinated axonal regions
float RM_my     = 10        // myelinated axon segments.
float ELEAK_sd  = -0.060    // soma & dend
float ELEAK_ax  = -0.060    // axon
float EREST_ACT = -0.060
