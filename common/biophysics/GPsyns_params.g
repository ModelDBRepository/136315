// ===============================================================
//            Synaptic conductances
// ===============================================================

// STN AMPA inputs
float G_AMPA        = 0.25e-9    
float tauRise_AMPA  = 0.001
float tauFall_AMPA  = 0.003

// Striatal inhibitory inputs
float G_GABA        = 0.25e-9
float tauRise_GABA  = 0.001
float tauFall_GABA  = 0.012

// Default input rates = 0
float STN_rate      = 0
float striatum_rate = 0

// Random seeds for timetables
float rseed_STN     = 78923456
float rseed_Str     = 78123456

// Reversal potentials
float E_AMPA        = 0
float E_GABA        = -0.080

