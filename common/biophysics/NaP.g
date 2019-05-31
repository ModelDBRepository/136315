//==================================================================
//                      Persistent Na channel
//    Based on Magistretti & Alonso (1999), JGP 114:491-509
//        and Magistretti & Alonso (2002), JGP 120: 855-873.
//    Created by J.R. Edgerton, 03/2004
//    Modified 10/2004 by JRE: add z-gate slow inactivation, improve 
//        model's y-gate intermediate inactivation.
//==================================================================

// --> kinetics for room temperature
float mpower_NaP    = 3
float Vhalfm_NaP    = -0.050
float Km_NaP        = 0.0057
float taummin_NaP   = 0.00003        // room temperature
float taummax_NaP   = 0.00043686    // room temperature
float V0taum_NaP    = -0.04264
float Ktaum1_NaP    = 0.0144
float Ktaum2_NaP    = -0.0144

float hpower_NaP    = 1
float hmin_NaP      = 0.154
float V0h_NaP       = -0.057
float Kh_NaP        = -0.004
float tauhmin_NaP   = 0.03            // room temp
float tauhmax_NaP   = 0.051            // room temp
float V0tauh_NaP    = -0.034
float Ktauh1_NaP    = 0.026
float Ktauh2_NaP    = -0.0319

// Couldn't get the same curve shape with the standard tau(V) equation.
float spower_NaP    = 1
float V0s_NaP       = -0.01
float Ks_NaP        = -0.0049
float Aalpha_NaP    = -2.88        // units of /volt/sec
float Balpha_NaP    = -0.049    // units of /sec
float Kalpha_NaP    = 0.00463    // units of volts
float Abeta_NaP     = 6.94        // units of /volt/sec
float Bbeta_NaP     = 0.447        // units of /sec
float Kbeta_NaP     = -0.00263    // units of volts

float dq10_NaP      = 3                // divide all tau values by this number
function make_Na_slow_GP
    if ({exists Na_slow_GP})
        return
    end
    create  tabchannel Na_slow_GP
    setfield Na_slow_GP Ek {ENa} Gbar {G_Na_slow_GP} Ik 0 Gk 0 \
    Xpower {mpower_NaP} Ypower {hpower_NaP} Zpower {spower_NaP}
    
//    ***    Activation & Deactivation (m-gate)
    float Km        = {Km_NaP}
    float Vhalfm    = {Vhalfm_NaP}
    float V0taum    = {V0taum_NaP}
    float taummax   = {taummax_NaP} / {dq10_NaP}
    float taummin   = {taummin_NaP} / {dq10_NaP}
    float Ktau1     = {Ktaum1_NaP}
    float Ktau2     = {Ktaum2_NaP}
    float minf, taum
    
    float V0m = {Vhalfm} + ({Km} * {log {(1 / {pow 0.5 {1/{mpower_NaP}}}) - 1}})
    echo "Na_slow V0m: " {V0m}

    call Na_slow_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        minf = 1 / (1 + {exp {({V0m} - x) / {Km} }})
        taum = {taummin} + (({taummax} - {taummin}) / ({exp { ({V0m} - x) / {Ktau1} } } + {exp {({V0m} - x) / {Ktau2} }}))
        setfield Na_slow_GP X_A->table[{i}] {taum}
        setfield Na_slow_GP X_B->table[{i}] {minf}
        x = x + dx
    end

    tweaktau Na_slow_GP X
    call Na_slow_GP TABFILL X 6000 0
    setfield Na_slow_GP X_A->calc_mode {NO_INTERP}
    setfield Na_slow_GP X_B->calc_mode {NO_INTERP}

//    ***    Fast / Intermediate Inactivation (h-gate)
    float hmin    = {hmin_NaP}
    float V0h     = {V0h_NaP}
    float Kh      = {Kh_NaP}
    float V0tauh  = {V0tauh_NaP}
    float Ktauh1  = {Ktauh1_NaP}
    float Ktauh2  = {Ktauh2_NaP}
    float tauhmin = {tauhmin_NaP} / {dq10_NaP}
    float tauhmax = {tauhmax_NaP} / {dq10_NaP}

    float tauh, hinf
    call Na_slow_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        
        hinf = {hmin} + ((1 - {hmin}) / (1 + {exp {({V0h} - x) / {Kh} }}))
        tauh = ({tauhmin} + (({tauhmax} - {tauhmin}) / ({exp {({V0tauh} - x) / {Ktauh1}}} + {exp {({V0tauh} - x) / {Ktauh2}}})))
            
        setfield Na_slow_GP Y_A->table[{i}] {tauh}
        setfield Na_slow_GP Y_B->table[{i}] {hinf}
        x = x + dx
    end
    
    tweaktau Na_slow_GP Y
    call Na_slow_GP TABFILL Y 6000 0
    setfield Na_slow_GP Y_A->calc_mode {NO_INTERP}
    setfield Na_slow_GP Y_B->calc_mode {NO_INTERP}

//    *** Slow Inactivation (s-gate)    
    float Ks       = {Ks_NaP}
    float V0s      = {V0s_NaP}
    float Aalpha   = {Aalpha_NaP}
    float Balpha   = {Balpha_NaP}
    float Kalpha   = {Kalpha_NaP}
    float Abeta    = {Abeta_NaP}
    float Bbeta    = {Bbeta_NaP}
    float Kbeta    = {Kbeta_NaP}

    float alphas, betas, taus, sinf
    call Na_slow_GP TABCREATE Z {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        alphas = (({Aalpha} * x) + {Balpha}) / (1 - {exp {((x + ({Balpha} / {Aalpha})) / {Kalpha})}})
        betas = (({Abeta} * x) + {Bbeta}) / (1 - {exp {((x + ({Bbeta} / {Abeta})) / {Kbeta})}})
        
        taus = 1 / ({dq10_NaP} * ({alphas} + {betas}))
        sinf = 1 / (1 + {exp {({V0s} - x) / {Ks} }})
        setfield Na_slow_GP Z_A->table[{i}] {taus}
        setfield Na_slow_GP Z_B->table[{i}] {sinf}
        x = x + dx
    end
    tweaktau Na_slow_GP Z
    call Na_slow_GP TABFILL Z 6000 0
    setfield Na_slow_GP Z_A->calc_mode {NO_INTERP}
    setfield Na_slow_GP Z_B->calc_mode {NO_INTERP}
    setfield Na_slow_GP Z_conc 0  // Z-gate voltage-dependent
end
