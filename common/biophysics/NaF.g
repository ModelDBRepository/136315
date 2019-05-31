//==================================================================
//                      Fast Na channel
//    Activation and fast inactivation made to replicate resurgent
//        sodium current from Raman & Bean as closely as possible.
//    Slow inactivation gate added by J. Edgerton, 2004.
//    Support for voltage-dependent Z-gate by Cengiz Gunay, 2004
//==================================================================
// --> kinetics for 32 degrees C:
float mpower_NaF    = 3
float Vhalfm_NaF    = -0.0324
float Km_NaF        = 0.005
float taummin_NaF   = 0.000028
float taummax_NaF   = 0.000028
float Ktaum1_NaF    = 1    // irrelevant because taumax == taumin
float Ktaum2_NaF    = 1 // irrelevant same as Ktaum1

float hpower_NaF    = 1
float V0h_NaF       = -0.048
float Kh_NaF        = -0.0028
float tauhmin_NaF   = 0.0003
float tauhmax_NaF   = 0.016    
float V0tauh_NaF    = -0.043
float Ktauh1_NaF    = 0.01
float Ktauh2_NaF    = -0.005

float spower_NaF    = 1
float mins_NaF      = 0.15
float V0s_NaF       = -0.040
float Ks_NaF        = -0.0054
float tausmin_NaF   = 0.01
float tausmax_NaF   = 1    
float Ktaus1_NaF    = 0.0183
float Ktaus2_NaF    = -0.010

float dq10_NaF      = 1

function make_Na_fast_GP
    if ({exists Na_fast_GP})
            return
    end
    create  tabchannel Na_fast_GP
    setfield Na_fast_GP Ek {ENa} Gbar {G_Na_fast_GP} Ik 0 Gk 0\
        Xpower {mpower_NaF} Ypower {hpower_NaF} Zpower {spower_NaF}
        
//    Activation & Deactivation
    float Vhalfm     = {Vhalfm_NaF}
    float Km         = {Km_NaF}
    float taummax    = {taummax_NaF} / {dq10_NaF}
    float taummin    = {taummin_NaF}    / {dq10_NaF}
    float Ktau1      = {Ktaum1_NaF}
    float Ktau2      = {Ktaum2_NaF}
    float V0m, minf, taum
    V0m = {Vhalfm} + ({Km} * {log {(1 / {pow 0.5 {1/{mpower_NaF}}}) - 1}})
    echo "Na_fast V0m: " {V0m}
    call Na_fast_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        minf = 1 / (1 + {exp {({V0m} - x) / {Km} }})
        taum = {taummin} + (({taummax} - {taummin}) / ({exp { ({V0m} - x) / {Ktau1} } } + {exp {({V0m} - x) / {Ktau2} }}))
        setfield Na_fast_GP X_A->table[{i}] {taum}
        setfield Na_fast_GP X_B->table[{i}] {minf}
        x = x + dx
    end
    tweaktau Na_fast_GP X
    call Na_fast_GP TABFILL X 6000 0
    setfield Na_fast_GP X_A->calc_mode {NO_INTERP}
    setfield Na_fast_GP X_B->calc_mode {NO_INTERP}

//    Fast Inactivation
    float V0h        = {V0h_NaF}
    float V0tauh     = {V0tauh_NaF}
    float Kh         = {Kh_NaF}
    float tauhmax    = {tauhmax_NaF} / {dq10_NaF} 
    float tauhmin    = {tauhmin_NaF} / {dq10_NaF}
    float Ktauh1     = {Ktauh1_NaF}
    float Ktauh2     = {Ktauh2_NaF}
    float hinf, tauh
    call Na_fast_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        hinf = 1 / (1 + {exp {({V0h} - x) / {Kh} }})
        tauh = {tauhmin} + (({tauhmax} - {tauhmin}) / ({exp { ({V0tauh} - x) / {Ktauh1} } } + {exp {({V0tauh} - x) / {Ktauh2} }}))
        setfield Na_fast_GP Y_A->table[{i}] {tauh}
        setfield Na_fast_GP Y_B->table[{i}] {hinf}
        x = x + dx
        end
    tweaktau Na_fast_GP Y
    call Na_fast_GP TABFILL Y 6000 0
    setfield Na_fast_GP Y_A->calc_mode {NO_INTERP}
    setfield Na_fast_GP Y_B->calc_mode {NO_INTERP}


//    Slow Inactivation
//    Equations & params from Spampanato et al, 2003, except that
//        tausmin added to prevent segmentation faults due to
//        excessively small time constants at voltage extremes.
    float V0s        = {V0s_NaF}
    float V0taus     = {V0s_NaF}
    float Ks         = {Ks_NaF}
    float mins       = {mins_NaF}
    float Ktaus1     = {Ktaus1_NaF}
    float Ktaus2     = {Ktaus2_NaF}
    float tausmax    = {tausmax_NaF} / {dq10_NaF}
    float tausmin    = {tausmin_NaF} / {dq10_NaF}
    float sinf, taus
    
    call Na_fast_GP TABCREATE Z {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        sinf = {mins} + ((1-{mins}) / (1 + {exp {({V0s} - x) / {Ks} }}))
        taus = tausmin + ({tausmax} - {tausmin}) / ({exp {({V0taus} - x) / {Ktaus1}}} + {exp {({V0taus} - x) / {Ktaus2}}}) 

        setfield Na_fast_GP Z_A->table[{i}] {taus}
        setfield Na_fast_GP Z_B->table[{i}] {sinf}
        x = x + dx
    end
    tweaktau Na_fast_GP Z
    call Na_fast_GP TABFILL Z 6000 0
    setfield Na_fast_GP Z_A->calc_mode {NO_INTERP}
    setfield Na_fast_GP Z_B->calc_mode {NO_INTERP}
    setfield Na_fast_GP Z_conc 0  // Z-gate voltage-dependent
end

