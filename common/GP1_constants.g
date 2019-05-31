// GP1_constants - set constant param values for GP simulations.

// File Paths
str morph_fname          = "../common/morphol/GP1.p"

str allcompsfilename     = "../common/morphol/GP1allcompnames.asc"
int ncomps               = 585   // total # compartments in model. 

str dendfilename         = "../common/morphol/GP1dendritenames.asc"
int num_dendcomps        = 511   // total # dendritic compartments

str chanscale_fname      = "../common/morphol/GP1_somaGeoDist_dendcompts.asc"

//simulation run defaults
int is_hsolved           = 1     // use hines solver
str cellpath             = "/GP"
float dt                 = 1e-5  // time step, in seconds
float PI                 = 3.14159

//Voltage-gated ion channel reversal potentials
float ENa                = 0.050
float ECa                = 0.130
float EK                 = -0.090
float Eh                 = -0.03
