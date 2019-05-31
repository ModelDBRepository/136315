//=============================================================================
//                                  KCNQ 
// Activation kinetics: Gamper, Stockand, Shapiro (2003). J Neurosci 23: 84-95.
// GV curve, deact kinetics: Prole & Marrion (2004). Biophys J. 86: 1454-69.
//=============================================================================

// --> kinetics for 32 degrees (Q10=3 adjusted from paper)

float npower_KCNQ     = 4
float Vhalfn_KCNQ     = -0.0285
float Kn_KCNQ         = 0.0195    // Yields actual K = 15 when power = 4.
float taunmin_KCNQ    = 0.0067
float taunmax_KCNQ    = 0.100
float Ktaun1_KCNQ     = 0.035
float Ktaun2_KCNQ     = -0.025

float dq10_KCNQ       = 1

function make_KCNQ_GP
    if (({exists KCNQ_GP}))
        return
    end
    create tabchannel KCNQ_GP
    setfield KCNQ_GP Ek {EK} Gbar {G_KCNQ_GP} Ik 0 Gk 0\
        Xpower {npower_KCNQ} Ypower 0 Zpower 0
    float Vhalfn  = {Vhalfn_KCNQ}    // True Vhalf for channel activation
    float Kn      = {Kn_KCNQ}
    float taunmax = {taunmax_KCNQ} / {dq10_KCNQ}
    float taunmin = {taunmin_KCNQ} / {dq10_KCNQ}
    float K1tau   = {Ktaun1_KCNQ}
    float K2tau   = {Ktaun2_KCNQ}
    float V0n, ninf, taun, alpha, beta
    V0n = {Vhalfn} + ({Kn} * {log {(1 / {pow 0.5 {1/{npower_KCNQ}}}) - 1}})
    echo "KCNQ V0n: " {V0n}
    //V0n is Vhalf for each individual n gate.
    call KCNQ_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        ninf = 1 / (1 + {exp { ({V0n} - x) / {Kn} } } )
        taun = {taunmin} + (({taunmax} - {taunmin}) / ({exp { ({V0n} - x) / {K1tau} } } + {exp { ({V0n} - x) / {K2tau} }}))
        setfield KCNQ_GP X_A->table[{i}] {taun}
        setfield KCNQ_GP X_B->table[{i}] {ninf}
        x = x + dx
    end
    tweaktau KCNQ_GP X
    call KCNQ_GP TABFILL X 6000 0
    setfield KCNQ_GP X_A->calc_mode {NO_INTERP}
    setfield KCNQ_GP X_B->calc_mode {NO_INTERP}
end
